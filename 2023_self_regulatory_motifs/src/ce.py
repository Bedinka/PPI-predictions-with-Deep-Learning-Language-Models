#!python

import os, sys, re

def cmp(a, b):
    return (a > b) - (a < b) 

def getCommandOutput(command, exedir = None):
    if exedir != None:
        os.chdir( exedir )
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        raise RuntimeError('%s failed w/ exit code %d' % (command, err))
    return data

class AlignedStructure:
    pass

def findMostFit( data ):
    output = []
    for line in data.split("\n"):
        #Alignment length = 55 Rmsd = 0.00A Z-Score = 5.3 Gaps = 0(0.0%) CPU = 158 ms. Sequence identities = 100.0%
        if line[:len("Alignment")] == "Alignment":
            fields = line.split()
            length = float( fields[3] )
            rmsd = float( fields[6][:-1] )
            zscore = float( fields[9] )
            print(length, rmsd, zscore)
            aAligned = AlignedStructure()
            aAligned.length = length
            aAligned.rmsd = rmsd
            aAligned.zscore = zscore
            output.append( aAligned )
        if line[:len("     X1 = ")] == "     X1 = ":
            aAligned.X = [ float(x) for x in re.findall( "\(\s*(\-*\d+\.\d+)\)", line ) ]
        if line[:len("     Y1 = ")] == "     Y1 = ":
            aAligned.Y = [ float(x) for x in re.findall( "\(\s*(\-*\d+\.\d+)\)", line ) ]
        if line[:len("     Z1 = ")] == "     Z1 = ":
            aAligned.Z = [ float(x) for x in re.findall( "\(\s*(\-*\d+\.\d+)\)", line ) ]
    output.sort(lambda x, y:cmp(x.zscore,y.zscore), reverse = True)  # 10 -> 9 -> 8 ...
    return output

def translatePDB( filename, aAligned ):
    # ATOM   2692  N   LYS G  12      55.763  32.551  37.146  1.00 55.62           N
    # ATOM   2693  CA  LYS G  12      56.776  32.217  36.152  1.00 55.72           C
    # ATOM   2694  C   LYS G  12      56.408  31.077  35.204  1.00 55.02           C
    # 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
    # 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
    # 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
    # X1 = ( 1.000000)*Xorig + ( 0.000000)*Yorig + ( 0.000000)*Zorig + (    0.000000)
    # Y1 = (-0.000000)*Xorig + ( 1.000000)*Yorig + ( 0.000000)*Zorig + (   -0.000000)
    # Z1 = (-0.000000)*Xorig + ( 0.000000)*Yorig + ( 1.000000)*Zorig + (   -0.000000)

    output = []
    f = open( filename )
    for line in f:
        if line[:len("ATOM")] == "ATOM":
            x0 = float( line[31:38] )
            y0 = float( line[38:46] )
            z0 = float( line[46:54] )
            #print line[:-1] + "##" 
            #print x0, y0, z0
            x1 = aAligned.X[0] * x0 + aAligned.X[1] * y0 + aAligned.X[2] * z0 + aAligned.X[3] 
            y1 = aAligned.Y[0] * x0 + aAligned.Y[1] * y0 + aAligned.Y[2] * z0 + aAligned.Y[3] 
            z1 = aAligned.Z[0] * x0 + aAligned.Z[1] * y0 + aAligned.Z[2] * z0 + aAligned.Z[3] 
            output.append( line[:30] + "%8.3f" % x1 + "%8.3f" % y1 + "%8.3f" % z1 + line[54:-1] )
            #print line[:30] + "%8.3f" % x1 + "%8.3f" % y1 + "%8.3f" % z1 + line[54:-1]
    f.close()
    return output

def inverseTranslatePDB( filename, aAligned ):
    # ATOM   2692  N   LYS G  12      55.763  32.551  37.146  1.00 55.62           N
    # ATOM   2693  CA  LYS G  12      56.776  32.217  36.152  1.00 55.72           C
    # ATOM   2694  C   LYS G  12      56.408  31.077  35.204  1.00 55.02           C
    # 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
    # 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
    # 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
    # X1 = ( 1.000000)*Xorig + ( 0.000000)*Yorig + ( 0.000000)*Zorig + (    0.000000)
    # Y1 = (-0.000000)*Xorig + ( 1.000000)*Yorig + ( 0.000000)*Zorig + (   -0.000000)
    # Z1 = (-0.000000)*Xorig + ( 0.000000)*Yorig + ( 1.000000)*Zorig + (   -0.000000)
    a = aAligned.X[0]
    b = aAligned.X[1]
    c = aAligned.X[2]
    d = aAligned.Y[0]
    e = aAligned.Y[1]
    f = aAligned.Y[2]
    g = aAligned.Z[0]
    h = aAligned.Z[1]
    k = aAligned.Z[2]

    l = aAligned.X[3]
    m = aAligned.Y[3]
    n = aAligned.Z[3]

    inv_detA = 1.0 / ( a * ( e*k -f*h ) + b * ( f*g - k*d ) + c * ( d*h - e*g ) )

    A = ( e*k - f*h )
    B = ( f*g - d*k )
    C = ( d*h - e*g )
    D = ( c*h - b*k )
    E = ( a*k - c*g )
    F = ( g*b - a*h )
    G = ( b*f - c*e )
    H = ( c*d - a*f )
    K = ( a*e - b*d )

    output = []
    f = open( filename )
    for line in f:
        if line[:len("ATOM")] == "ATOM":
            x1 = float( line[31:38] ) - l
            y1 = float( line[38:46] ) - m
            z1 = float( line[46:54] ) - n
            #print line[:-1] + "##" 
            #print x0, y0, z0
            x0 = inv_detA * ( A * x1 + D * y1 + G * z1 )
            y0 = inv_detA * ( B * x1 + E * y1 + H * z1 )
            z0 = inv_detA * ( C * x1 + F * y1 + K * z1 )
            output.append( line[:30] + "%8.3f" % x0 + "%8.3f" % y0 + "%8.3f" % z0 + line[54:-1] )
            #print line[:30] + "%8.3f" % x1 + "%8.3f" % y1 + "%8.3f" % z1 + line[54:-1]
    f.close()
    return output

if __name__ == "__main__":
    if len( sys.argv ) < 3:
        print("python ce.py ref query")
        sys.exit(0)

    pdb1 = sys.argv[1]
    pdb2 = sys.argv[2]
    if len( sys.argv ) == 4:
        outfile = os.path.abspath( sys.argv[3] )
    else:
        outfile = os.path.abspath( "newPDB.pdb" )
        inv_outfile = os.path.abspath( "inv_newPDB.pdb" )

    pdb1 = os.path.abspath( pdb1 )
    pdb2 = os.path.abspath( pdb2 )

    data = getCommandOutput( "/home/nisusers/jyang/bin/ce_java/runCE.sh -file1 %s -file2 %s -printCE" % ( pdb1, pdb2 ), "/home/nisusers/jyang/bin/ce_java" )

    aligned_list = findMostFit( data )

    print(aligned_list[0])
    print(aligned_list[0].X)
    print(aligned_list[0].Y)
    print(aligned_list[0].Z)

    newPDB = translatePDB( pdb2, aligned_list[0] )
    f = open( outfile, "w" )
    for line in newPDB:
        f.write( line )
        f.write( "\n" )
    f.close()

    inv_newPDB = inverseTranslatePDB( pdb1, aligned_list[0] )
    f = open( inv_outfile, "w" )
    for line in inv_newPDB:
        f.write( line )
        f.write( "\n" )
    f.close()


