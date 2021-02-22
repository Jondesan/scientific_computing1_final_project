#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 23:47:39 2020

@author: joonahuh
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)


#   TO DO:
#   Translate comments and variable names from FIN to ENG



#   Datan sisäänluku
'''
    5 saraketta:
    Aika (tunteina)
    Pienhiukkaspitoisuus (hiukkasia/cm^3)
    Pienhiukkasten häviönopeus (hiukkasia/cm^3/s)
    Rikkihapon H2SO4 pitoisuus (molekyyliä/cm^3)
    Orgaanisten yhdisteiden C10 pitoisuus (molekyyliä/cm^3)
'''
inptData = np.loadtxt('measurement_data.dat')
steadyDat = np.loadtxt('experiment_steady.dat')

#   Funktio luo datasta liukuvan keskiarvon
def liukuva_keskiarvo(data, n):
    
    datOut = np.array([])
    
    #   Käydään koko data läpi ja muutetaan liukuvaan keskiarvoon otettavien lukujen
    #   määrää riippuen ollaanko alussa, lopussa vai keskellä datalistaa
    for i in range( len(data) ):
        if i < (n-1)/2:
            sum = 0
            for j in range( int( ( n + 1 ) / 2 ) + i ):
                sum += data[j]
            avg = sum / ( ( n + 1 ) / 2 + i )
            datOut = np.append(datOut, avg)
        elif i + (n-1) / 2 > len(data) - 1:
            sum = 0
            
            for j in range( int( ( n - 1 ) / 2 - i + len(data) ) ):
                sum += data[ - ( j + 1 ) ]
            avg = sum/( (n-1)/2 - i + len(data) )
            datOut = np.append(datOut, avg)
            
        else:
            sum = 0
            for j in range(n):
                sum += data[ i + j - int( ( n - 1 ) / 2 ) ]
            avg = sum/n
            datOut = np.append(datOut, avg)
    #   Palautetaan muokattu datataulukko
    return datOut


#   Funktio palauttaa sisään syötettyjen datojen approksimoidun derivaatan
def derivaatta(xDat, yDat):
    
    #   Jos listat erisuuret, palauta NaN
    if len(xDat) != len(yDat):
        return np.nan
    
    derivative = np.array([])
    
    #   Käydään listat läpi ja lasketaan derivaatan approksimaatiot
    #   Alussa ja lopussa lisätään NaN.
    for i in range(len(xDat)):
        if i == 0 or i == (len(xDat) - 1):
            derivative = np.append(derivative, np.nan)
        else:
            value = (yDat[i+1] - yDat[i-1]) / (xDat[i+1] - xDat[i-1])
            derivative = np.append(derivative, value)
    
    return derivative

# Funktio palauttaa datasarjan steady state -jakson keskiarvot listassa
def steady_data(dataList1):
    
    time = inptData[:,0]
    data = np.array([])
    steadyUpperLim = steadyDat[:,1]
    steadyLowerLim = steadyDat[:,0]
    
    for i in range( len( steadyDat[:,0] ) ):        
        
        temp = np.array([])
        
        for j in range( len(time) ):
            if time[j] >= steadyLowerLim[i] and time[j] <= steadyUpperLim[i]:
                temp = np.append(temp, dataList1[j])
        data = np.append(data, np.nanmean(temp))
    return data


#   Tehdään järkevät muuttujat datan liukuvalle keskiarvolle
particleDat = liukuva_keskiarvo(inptData[:,1],5)
sulfuricDat = liukuva_keskiarvo(inptData[:,3],5)
orgDat = liukuva_keskiarvo(inptData[:,4],5)
jDat = derivaatta(inptData[:,0]*(60**2), particleDat) + liukuva_keskiarvo(inptData[:,2],5)


#   Tehdään tarvittavat kuvaajat subplotteina
fig1, ax1 = plt.subplots(3,1, figsize=(12,6))
fig2, ax2 = plt.subplots(figsize=(8,4))
fig3, ax3 = plt.subplots(figsize=(12,8))
fig4, ax4 = plt.subplots(figsize=(12,8))
fig5, ax5 = plt.subplots(figsize=(12,8))


ax1[0].semilogy(inptData[:,0],particleDat)
ax1[1].semilogy(inptData[:,0],sulfuricDat)
ax1[2].semilogy(inptData[:,0],orgDat)


for i in range(3):
    ax1[i].set_xlim(xmin=0, xmax=120)
    ax1[i].set_xlabel('Aika (tunteina)')
ax2.set_xlim(xmin=75, xmax=78)
ax2.set_ylim(ymin=4,ymax=170)
ax2.set_xlabel('Aika (tunteina)')



ax1[0].set_ylabel('Keskiarvoistettu\n pienhiukkaspitoisuus $1/{\mathrm{cm}^3}$')
ax1[1].set_ylabel('Keskiarvoistettu\n rikkihappopitoisuus $1/{\mathrm{cm}^3}$')
ax1[2].set_ylabel('Keskiarvoistettu orgaanisten\n yhdisteiden pitoisuus $1/{\mathrm{cm}^3}$')

ax2.semilogy(inptData[:,0],
             inptData[:,1],
             label='Pienhiukkaspitoisuus')

ax2.semilogy(inptData[:,0],
             particleDat,
             label='Keskiarvoistettu pienhiukkaspitoisuus')

ax2.set_ylabel('Pienhiukkaspitoisuus $1/{\mathrm{cm}^3}$')



ax3.loglog(sulfuricDat,
           jDat,
           linestyle='none',
           marker='.',
           markersize=8,
           color='blue')
ax3.set_xlabel('Rikkihappopitoisuus $\mathrm{H_2SO_4}$ [$1/\mathrm{cm^3}$]')
ax3.set_ylabel('Pienhiukkasten muodostumisnopeus $J$ [$1/\mathrm{cm^3s}$]')


ax4.loglog(orgDat,
           jDat,
           linestyle='none',
           marker='.',
           markersize=8,
           color='red')
ax4.set_xlabel('Orgaanisten yhdisteiden $\mathrm{C_{10}}$ pitoisuus [$1/\mathrm{cm^3}$]')
ax4.set_ylabel('Pienhiukkasten muodostumisnopeus $J$ [$1/\mathrm{cm^3s}$]')


ax5.loglog(orgDat * sulfuricDat,
           jDat,
           linestyle='none',
           marker='.',
           markersize=8,
           color='green')
ax5.set_xlabel('Orgaanisten yhdisteiden $\mathrm{C_{10}}$ ja rikkihapon $\mathrm{H_2SO_4}$'
               'pitoisuuksien tulo [$1/\mathrm{cm^3}$]')
ax5.set_ylabel('Pienhiukkasten muodostumisnopeus $J$ [$1/\mathrm{cm^3s}$]')




ax2.legend(loc='lower right')

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()
fig1.savefig('kuva1.png',
            dpi=300)
fig2.savefig('kuva2.png',
             dpi=300)
fig3.savefig('kuva3.png',
             dpi=300)
fig4.savefig('kuva4.png',
             dpi=300)
fig5.savefig('kuva5.png',
             dpi=300)


#   Tehdään datalle steady state taulukot
steadyStateSulfuric = steady_data(sulfuricDat)
steadyStateOrganic = steady_data(orgDat)
steadyJDat = steady_data( jDat )

fig6, ax6 = plt.subplots(figsize=(8,4))


ax6.loglog(steadyStateSulfuric * steadyStateOrganic,
           steadyJDat,
           linestyle='none',
           marker='.',
           markersize=8,
           color='black',
           label='Steady-state datapisteet')

x = np.log10(steadyStateSulfuric * steadyStateOrganic)
y = np.log10(steadyJDat)
params = np.polyfit(x,y,1)

xx = np.logspace(np.min(x),
                 np.max(x),
                 50)
yy = 10**params[1] * xx**params[0]

ax6.plot(xx,
         yy,
         color='red',
         linewidth=4,
         label='Steady-state sovitus')

ax6.set_xlabel('Rikkihapon ja orgaanisten yhdisteiden pitoisuuksien steady\n'
               'state -tilojen tulon 10-kantainen logaritmi [$\log_{10}({1/\mathrm{cm}^3})$]')
ax6.set_ylabel('Hiukkasmuodostuksen $J$ steady state -tilojen keskiarvon\n'
               '10-kantainen logaritmi [$\log_{10}(1/\mathrm{cm^3s})$]')


ax6.legend(loc='upper left')

fig6.tight_layout()
fig6.savefig('kuva6.png',
             dpi=300)







