from scitools.std import *
import sys

def plot_orbit():
    infile1 = open('positions.dat', 'r')

    # number of objects
    n = int(infile1.readline())

    # read fil and append to lists
    # x30 is the third x-value of object zero:
    # positions.dat is made in the following way: x00 y00 z00 x01 y01 z01 x02 y02 z02 ... x0n-1 y0n-1 z0n-1
    #                                             x10 y10 z10 x11 y11 z11 x12 y12 z12 ... x1n-1 y1n-1 z1n-1
    #                                             ..   ..  ..  ..  ..  .. ...  ..     ..
    #                                             xN0 yN0 zN0 xN1 yN1 zN1 xN2 yN2 zN2 ... xNn-1 yNn-1 zNn-1
    # I will make three lists x, y and z: 
    # x = [[x00, x01, x02, ... , x0n-1], [x10, x11, x12, ... , x1n-1], ... , [xN0, xN1, xN2, ... , xNn-1]]
    # ... and similar for y and z 
                     
    x = []
    y = []
    z = []
    for line in infile1:
        words = line.split()
        xx = []
        yy = []
        zz = []
        for i in range(n):
            xx.append(float(words[3*i]))
            yy.append(float(words[3*i+1]))
            zz.append(float(words[3*i+2]))
        x.append(xx)
        y.append(yy)
        z.append(zz)

        infile1.close()

    x = array(x) 	# [[x00 x01 x02], [x10 x11 x12] ...
    y = array(y)
    z = array(z)

    # the x-values of the Sun is then the first column of x, which is indexed as x[::,0] and so on
    for i in range(n):
        plot(x[::,i], y[::,i])
        hold('on')
    hold('off')

    if n == 2:
        legend('Sun', 'Earth')
    elif n == 3:
        legend('Sun', 'Earth', 'Jupiter')
    elif n == 10:
        legend('Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto')

    xlabel('x position [AU]')
    ylabel('y position [AU]')
    axis('equal')

    raw_input('Press Return key to quit: ')


def plot_energymom():
    infile2 = open('energymom.dat', 'r')

    t_max = int(infile2.readline())

    K = []
    U = []
    L = []
    for line in infile2:
        words = line.split()
        K.append(float(words[0]))
        U.append(float(words[1]))
        L.append(float(words[2]))

    infile2.close()

    K = array(K)
    U = array(U)
    L = array(L)

    subplot(2,2,1)
    plot(array(range(len(K)))/t_max, K, 'b-')
    legend('K')
    xlabel('Time [yr]')
    ylabel(r'Energy [$\frac{AU^2}{yr^2}$]')
    subplot(2,2,2)
    plot(array(range(len(K)))/t_max, U, 'g-')
    xlabel('Time [yr]')
    ylabel(r'Energy [$\frac{AU^2}{yr^2}$]')
    legend('U')
    subplot(2,2,3)
    plot(array(range(len(K)))/t_max, L, 'y-')
    xlabel('Time [yr]')
    ylabel(r'Ang. mom. [$\frac{AU}{yr^2}$]')
    legend('Lz')
    subplot(2,2,4)
    plot(array(range(len(K)))/t_max, K+U, 'r-')
    xlabel('Time [yr]')
    ylabel(r'Energy [$\frac{AU^2}{yr^2}$]')
    legend('K+U')
    raw_input('Press Return key to quit: ')


def extract_radius():
    infile3 = open('radius.dat', 'r')

    t_max = float(infile3.readline())
    r = float(infile3.readline())   
    
    return r, t_max
        

def plot_radius(method):
    import os
    no_of_values = 20
    N_values = linspace(500, 10000, no_of_values)

    final_radius = zeros(no_of_values)
    dt = zeros(no_of_values)
    for i in xrange(no_of_values):
        N = N_values[i]
        os.system('c++ main.cpp -larmadillo -llapack -lblas planet.h planet.cpp constants.h constants.cpp')
        os.system('./a.out %d %d' % (N, method))
        r, t_max = extract_radius()
        final_radius[i] = r
        dt[i] = t_max/(N+1)

    plot(dt, final_radius)    
    raw_input('Press Return key to quit: ')


### main ###
#plot_orbit()
#plot_energymom()'

method = int(sys.argv[1])
plot_radius(method)



