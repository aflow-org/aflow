#!/usr/bin/python3

Name = 'aurostd'
H = '.h'
C = '.cpp'

def convolve( type, file, includes ):
    fh = open( file, 'r' )
    for line in fh:
        if( line.startswith( '#include' ) and '"' in line):
            newfile = line.split('"')[1].strip()
            if( type is C and newfile == Name+H and not includes ):
                print(line, end='', file=oh)
            else:
                if( newfile.startswith('..') ):
                    print('#include "' + newfile[3:] + '"', file=oh)
                else:
                    if( newfile not in includes ):
                        includes.append( newfile )
                        convolve( type, newfile, includes )
        else:
            print(line, end='', file=oh)

oh = open('../'+Name+H, 'w')
convolve( H, Name+H, [] )
oh.close()

oh = open('../'+Name+C, 'w')
convolve( C, Name+C, [] )
oh.close()

