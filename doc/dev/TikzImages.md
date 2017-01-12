@page TikzImages Compiling the TikZ Images

If you want to update the TikZ images (mostly in the InterpolationTypes), you
need recompile the `*.tikz` files and convert them to PNG. You'll find them in
`doc/images`. To compile to PDF use

    pdflatex Tetrahedron3D.tikz

and then convert them to PNG using

    convert -density 200 Tetrahedron3D.pdf Tetrahedron3D.png

Thank you for flying with LaTeXTikZImageMagick Air.
