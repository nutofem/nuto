#pragma once
#include <boost/iterator/iterator_facade.hpp>

namespace NuTo
{
//! Iterator that points to values of TEntry produced `lazily` by TEntryGenerator. It implements the
//! boost::iterator facade http://www.boost.org/doc/libs/1_66_0/libs/iterator/doc/iterator_facade.html
//! @tparam TEntryGenerator generator type that provides the methods
//!    begin() / end()
//!    IsValid()
//!    Next()
//!    Get()
//! @tparam TEntry entry type
template <typename TEntryGenerator, typename TEntry>
class EntryIterator : public boost::iterator_facade<EntryIterator<TEntryGenerator, TEntry>, TEntry,
                                                    boost::single_pass_traversal_tag>
{
public:
    //! empty ctor indicates that the generator TEntryGenerator reached its end
    EntryIterator()
        : g(nullptr)
    {
    }

    //! ctor with a entry generator
    EntryIterator(TEntryGenerator* g)
        : g(g)
    {
    }

private:
    friend class boost::iterator_core_access;
    TEntry& dereference() const
    {
        if (!g || !g->IsValid())
            throw std::runtime_error("...");
        return g->Get();
    }
    bool equal(EntryIterator const& l) const
    {
        return g == l.g || ((!g || !g->IsValid()) && (!l.g || !l.g->IsValid()));
    }
    void increment()
    {
        if (g)
            g->Next();
    }

    TEntryGenerator* g;
};

//! Interface for a generator that generates an TEntry when calling Get()
//! @tparam TEntry entry type
template <typename TEntry>
class EntryGenerator
{
public:
    using TEntryIterator = EntryIterator<EntryGenerator, TEntry>;

    virtual ~EntryGenerator() = default;

    TEntryIterator begin()
    {
        return TEntryIterator(this);
    }
    TEntryIterator end()
    {
        return TEntryIterator();
    }

    //! Indicates the end of the generator
    //! @return false if there are no more TEntryIterators to generate
    virtual bool IsValid() const = 0;

    //! Modifies the internal generator state to move to the next entry
    virtual void Next() = 0;

    //! Generate a new entry
    //! @return iterator to the new entry
    virtual TEntry& Get() = 0;
};
} /* NuTo */
