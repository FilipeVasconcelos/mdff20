#!/usr/bin/python3
# ===============================================================
# version mdff    : 30 nov 2016 
# version mdff20  : avril 2020
# author          : fmv
# description     : This script generate gnuplot plots from input OSZIFF file
# ===============================================================

import matplotlib.pyplot as plt
import numpy as np
import datetime
import argparse
import sys


name_quant=['step','Time','Etot','Ekin','Utot','U_vdw','U_coul','Temp','Press','Pvir_vdw','Pvir_coul','Volume','Htot']
units={}
units['Time']='[ps]'
units['Etot']='[eV]'
units['Ekin']='[eV]'
units['Utot']='[eV]'
units['U_vdw']='[eV]'
units['U_coul']='[eV]'
units['Temp']='[K]'
units['Press']='[GPa]'
units['Volume']='[${\AA}^3$]'
units['Htot']='[eV]'

def read_OSZIFF(filename):

    alldata=[];step=[];time=[];etot=[];ekin=[];utot=[];uvdw=[];ucou=[]
    temp=[];pres=[];pvir_vdw=[];pvir_coul=[];volu=[];htot=[]

    k=0
    f=open(filename)
    for line in f :
        l = line.split()
        if len(l)==7 and  l[0][0] != "-" and l[0][0] !="s" :
            k+=1
            if k%2 != 0:
                step.append(int(l[0]))
                time.append(float(l[1]))
                etot.append(float(l[2]))
                ekin.append(float(l[3]))
                utot.append(float(l[4]))
                uvdw.append(float(l[5]))
                ucou.append(float(l[6]))
            if k%2 == 0:
                temp.append(float(l[1]))
                pres.append(float(l[2]))
                pvir_vdw.append(float(l[3]))
                pvir_coul.append(float(l[4]))
                volu.append(float(l[5]))
                htot.append(float(l[6]))
    if k%2 != 0:
        step.pop()
        time.pop()
        etot.pop()
        ekin.pop()
        utot.pop()
        uvdw.pop()
        ucou.pop()
    f.close()

    alldata.append(step)
    alldata.append(time)
    alldata.append(etot)
    alldata.append(ekin)
    alldata.append(utot)
    alldata.append(uvdw)
    alldata.append(ucou)
    alldata.append(temp)
    alldata.append(pres)
    alldata.append(pvir_vdw)
    alldata.append(pvir_coul)
    alldata.append(volu)
    alldata.append(htot)

    return alldata
    
def averaging(name_quant,alldata,last_points):

    average=[]
    average.append(None)#step
    average.append(None)#time
    for i,l in enumerate(alldata):
        if name_quant[i] != "step" and name_quant[i] != "Time":
            average.append(np.mean(l[-last_points:]))
            print ('{0:<11} {1:<10} {2:15.8e} {3:^10} {4:15.8e} '.format("<"+name_quant[i]+">","=",np.mean(l[-last_points:]),"std.",np.std(l[-last_points:])))

    return average

def plot_quant2(alldata,name_quant,average,q,title,l):

    t = alldata[1][-l:] 
    s = alldata[q[0]][-l:] 
    x1= [average[q[0]]]*len(t)
    plt.plot(t, s,'b')
    plt.plot(t, x1, '--' , color='b') 

    s = alldata[q[1]][-l:] 
    x1= [average[q[1]]]*len(t)
    plt.plot(t, s , 'g')
    plt.plot(t, x1, '--', color='g') 

    plt.xlabel('Time'+' '+units['Time'])
    plt.ylabel(name_quant[q[0]]+' '+units[name_quant[q[0]]])
    plt.title(title)
    plt.show()

def plot_quant(alldata,name_quant,average,q,title,l):

    t = alldata[1][-l:] 
    s = alldata[q][-l:] 

    x1= [average[q]]*len(t)
    plt.plot(t, s,'b')
    plt.plot(t, x1, '--' , color='b') 

    plt.xlabel('Time'+' '+units['Time'])
    plt.ylabel(name_quant[q]+' '+units[name_quant[q]])
    plt.title(title)
    plt.show()


def main_parser():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i','--input', help='Input file name',required=True,dest="input_file")
    parser.add_argument("-n", "--no_plot",dest="plot_flag",action="store_false",default=True,
            help="show plot of average and instantaneous thermodynamic quantities along OSZIFF")
    parser.add_argument("-l", "--last",dest="last_points",default=None,
            help="averaging on last <l> points")
    parser.add_argument("-v", "--volume",dest="lvolume",action="store_true",default=False,
            help="plot volume data")
    parser.add_argument("-p", "--pressure",dest="lpress",action="store_true",default=False,
            help="plot pressure data")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    separator=60*"="
    now = datetime.datetime.now()

    print( separator)
    print( now.strftime("%Y-%m-%d %H:%M"))
    print( "author : filipe.manuel.vasconcelos@gmail.com")
    print( separator)
    print( "Running poszi ...")
    print( "This script generate gnuplot plots from input OSZIFF file")

    args = main_parser()

    input_file=args.input_file
    lpress=args.lpress
    lvolume=args.lvolume

    alldata=read_OSZIFF(input_file)
    

    if args.last_points == None :
        last_points = len(alldata[0])
    else:
        last_points = int( args.last_points ) 

    print( len(alldata[0])," points in input file")
    print( "averaging and plot on last",last_points,"points")

    average=averaging(name_quant,alldata,last_points)
   
    if args.plot_flag :
        plot_quant2(alldata,name_quant,average,[2,12],title='Total Energy MD',l=last_points)
        plot_quant(alldata,name_quant,average,4,title='Potential Energy MD',l=last_points)
        plot_quant(alldata,name_quant,average,7,title='Temperature MD',l=last_points)
        if (lpress): plot_quant(alldata,name_quant,average,8,title='Pressure',l=last_points)
        if (lvolume): plot_quant(alldata,name_quant,average,11,title='Volume MD',l=last_points)
