function(warning)
    string(ASCII 27 Esc)
    set(ColourReset "${Esc}[m")
    set(Red         "${Esc}[31m")
    message(STATUS "${Red}${ARGV}${ColourReset}")
endfunction()

