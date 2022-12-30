# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 01:20:17 2022

@author: sameer_poudel
"""
from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

import numpy as np
import matplotlib.pyplot as plt
from functions_v5 import get_R_vs_Theta_Rb_values, get_UniaxYieldStress_vs_Theta_values, get_stress_points_in_Yieldlocus

class widgets(Frame):
    
    def __init__(self, window, picks=[]):
       
      Frame.__init__(self, window)
      
      self.checks = []
      self.curves = []
      self.curvenames = ['C1','C2','C3','C4','C5']
      lf = LabelFrame(window, text='Parameters')
      lf.grid(column=0, row=0)
      #lf.grid_propagate(True)
      
      for j in range(len(self.curvenames)): 
          
          label = Label(lf, text=self.curvenames[j], width=6)
          label.grid(row=j+1, column=0)
          
          self.vars = []
          for i, pick in enumerate(picks):
              
              if j == 0:
                  label = Label(lf, text=pick, width=6)
                  label.grid(row=j, column=i+1)
              
              if pick != 'select':       
                  var = IntVar()                            
                  if pick == 'e':
                    var.set(8)
                  else:
                    var.set(1)
                            
                  ent = Entry(lf, textvariable=var, width=6)
                  ent.grid(row=j+1, column=i+1)   
                  self.vars.append(ent)
                  
          chk = IntVar()
          cb = Checkbutton(lf, variable=chk)
          cb.grid(row=j+1,column=len(picks))
          
          self.curves.append(self.vars)
          self.checks.append(chk)

      Button(lf, text='Get the plot!', command=allstates).grid(row=(len(self.curvenames)+1),column=1)
      Button(lf, text='Normalize the plot!', command=normalize).grid(row=(len(self.curvenames)+1),column=3)   
      
    def state(self): 
      result_cvs = []
      result_names =  []
      result_chks = list(map((lambda chk: chk.get()), self.checks))
      for k, check in enumerate(result_chks):
          if check ==1:
              cvs = list(map((lambda ent: float(ent.get())), self.curves[k]))
              result_cvs.append(cvs)
              result_names.append(self.curvenames[k])
      return result_cvs, result_names
  
    
class display(Frame):
    
    def __init__(self, window, values, curve_names, picks=[]):
       
      Frame.__init__(self, window)
      
      lf = LabelFrame(window, text='Values')
      lf.grid(column=0, row=15)
      #lf.grid_propagate(True)
      
      for j, value in enumerate(values):
          
          label = Label(lf, text=curve_names[j], width=6)
          label.grid(row=j+1, column=0)
          
          #self.vars = []
          
          for i, pick in enumerate(picks):
              
              if j == 0:
                  label = Label(lf, text=pick, width=6)
                  label.grid(row=j, column=i+1)
              
              var =  IntVar(value=value[i])
              ent = Entry(lf, textvariable=var, width=6)
              ent.grid(row=j+1, column=i+1)   

                   
def allstates():
     
    global count, canvas, fig, ax1, ax2, ax3, toolbar
    
    count +=1 
    
    if count != 1:
        ax1.cla()
        ax2.cla()
        ax3.cla()
        
    curves, curve_names = lng.state()  
    print(curves)
    
    deg = np.array([list(range(0,91))])
    #ax2 = ax1.twinx()
    values = []
    
    for l, curve in enumerate(curves):
        
        R, Rb = get_R_vs_Theta_Rb_values(curve[:-1], curve[-1])
        stress = get_UniaxYieldStress_vs_Theta_values(curve[:-1], curve[-1], Y_0=1.0, NormalizeAlpha=False)
        #ax1.scatter(deg, R, marker='.')       
        ax1.scatter(deg, R, marker='d')
        
        ax2.scatter(deg,stress, marker='.', linewidths=0.0)
        #ax2.plot(deg,stress)
        #ax2.plot(deg,stress)
        yl = get_stress_points_in_Yieldlocus(curve[:-1],curve[-1],shear=0)
        ax3.plot(yl[0][0][0],yl[0][1][0], label = curve_names[l])
        
        ax3.scatter(yl[0][2],0,s=70,marker='^',color='#FF3339',label='Sigma_00',zorder=2)
        ax3.scatter(0,yl[0][3],s=70,marker='^',color='#F615DE',label='Sigma_90',zorder=2)
        ax3.scatter(yl[0][4],yl[0][4],s=70,marker='^',color='#1521F6',label='Sigma_45',zorder=2)
        ax3.scatter(yl[0][5],yl[0][5],s=70,marker='^',color='#15F669',label='Sigma_biax',zorder=2) 
        ax3.scatter(- yl[0][6],yl[0][6],s=70,marker='^',color='#F3F615',label='SigmaShear1',zorder=2)
        ax3.scatter(yl[0][7], - yl[0][7],s=70,marker='^',color='#F6AE15',label='SigmaShear2',zorder=2)
        ax3.scatter(yl[0][8],yl[0][9],s=70,marker='^',color='#070707',label='SigmaPlaneStrain_1',zorder=2)
        ax3.scatter(yl[0][10],yl[0][11],s=70,marker='^',color='#60A380',label='SigmaPlaneStrain_2',zorder=2) 
        
        #print(yl)
        
        v = [ round(float(yl[0][2]),4), round(float(yl[0][4]),4), round(float(yl[0][3]),4), round(float(yl[0][5]),4), 
                   round(float(yl[0][8]),4), round(float(yl[0][9]),4), round(float(yl[0][10]),4),  round(float(yl[0][11]),4),
                  round(float(yl[0][6]),4), round(float(yl[0][7]),4),
                  round(float(R[0][0]),4),  round(float(R[0][45]),4), round(float(R[0][90]),4), round(float(Rb),4)]
        
        values.append(v)
        
    display(frame1, values, curve_names, ['S00', 'S45', 'S90','Sbiax', 'PS1x', 'PS1y', 'PS2x', 'PS2y', 'S1', 'S2', 'R00',  'R45', 'R90', 'Rbiax'])    

    ax1.set_ylabel('R-values \u2666 ')
    ax1.set_xlabel('angle wrt RD')  
    ax1.set_title('Parameter Sensitivity Plot for Lankford Anistropy Parameters')      
    ax2.set_ylabel('Stress-values \u2022')
    ax3.set_title('Parameter Sensitivity Plot in Principal Stress Space')
    
    ax3.set_xlabel('sigma_xx')
    ax3.set_ylabel('sigma_yy')
    ax3.grid()
    

    #ax3.legend(loc='center')
    ax3.legend(bbox_to_anchor=(1.5, 1.0), loc='upper right')
        
    canvas.draw()
    canvas.get_tk_widget().grid(row = 0, column=15)
    #canvas.get_tk_widget().grid(row =5, column=0)
    
    
def normalize():
     
    global count, canvas, fig, ax1, ax2, ax3, toolbar
    
    count +=1 
    
    if count != 1:
        ax1.cla()
        ax2.cla()
        ax3.cla()
        
    curves, curve_names = lng.state()  
    #print(curves)
    
    deg = np.array([list(range(0,91))])
    #ax2 = ax1.twinx()
    
    values = []
    for l, curve in enumerate(curves):
        
        alpha = curve[:-1]
        print(alpha)
        a = curve[-1]
        alpha1 = ( ( (2 - abs((2*alpha[2]-2*alpha[3])/3.0)**a - abs((4*alpha[4]-alpha[5])/3)**a)**(1/a) )  * 3.0 - alpha[1] ) / 2.0
        
        if np.iscomplex(alpha1)==True:
            flag = 'complex'           
        else:
            flag = 'real'
            
        alpha[0] = alpha1.real
        print(alpha)
            
        R,Rb = get_R_vs_Theta_Rb_values(alpha, a)
        stress = get_UniaxYieldStress_vs_Theta_values(curve[:-1], curve[-1], Y_0=1.0, NormalizeAlpha=True)
        #print(stress)
        ax1.scatter(deg, R, marker='.')
        ax2.scatter(deg,stress,marker='d')
        #ax2.plot(deg,stress)
    
        yl = get_stress_points_in_Yieldlocus(alpha,a,shear=0)
        ax3.plot(yl[0][0][0],yl[0][1][0], label = f'{curve_names[l]}, a1 = {alpha[0]:.2f},{flag}' )
        
        ax3.scatter(yl[0][2],0,s=70,marker='^',color='#FF3339',label='Sigma_00',zorder=2)
        ax3.scatter(0,yl[0][3],s=70,marker='^',color='#F615DE',label='Sigma_90',zorder=2)
        ax3.scatter(yl[0][4],yl[0][4],s=70,marker='^',color='#1521F6',label='Sigma_45',zorder=2)
        ax3.scatter(yl[0][5],yl[0][5],s=70,marker='^',color='#15F669',label='Sigma_biax',zorder=2) 
        ax3.scatter(- yl[0][6],yl[0][6],s=70,marker='^',color='#F3F615',label='SigmaShear1',zorder=2)
        ax3.scatter(yl[0][7], - yl[0][7],s=70,marker='^',color='#F6AE15',label='SigmaShear2',zorder=2)
        ax3.scatter(yl[0][8],yl[0][9],s=70,marker='^',color='#070707',label='SigmaPlaneStrain_1',zorder=2)
        ax3.scatter(yl[0][10],yl[0][11],s=70,marker='^',color='#60A380',label='SigmaPlaneStrain_2',zorder=2) 
        
        
        v = [ round(float(yl[0][2]),4), round(float(yl[0][4]),4), round(float(yl[0][3]),4), round(float(yl[0][5]),4), 
                   round(float(yl[0][8]),4), round(float(yl[0][9]),4), round(float(yl[0][10]),4),  round(float(yl[0][11]),4),
                  round(float(yl[0][6]),4), round(float(yl[0][7]),4),
                  round(float(R[0][0]),4),  round(float(R[0][45]),4), round(float(R[0][90]),4), round(float(Rb),4)]
        
        values.append(v)
        
    display(frame1, values, curve_names, ['S00', 'S45', 'S90','Sbiax', 'PS1x', 'PS1y', 'PS2x', 'PS2y', 'S1', 'S2', 'R00',  'R45', 'R90', 'Rbiax'])    

    ax1.set_ylabel('R-values \u2666 ')
    ax1.set_xlabel('angle wrt RD')  
    ax1.set_title('Parameter Sensitivity Plot for Lankford Anistropy Parameters')      
    ax2.set_ylabel('Stress-values \u2022')
    ax3.set_title('Parameter Sensitivity Plot in Principal Stress Space')
    
    ax3.set_xlabel('sigma_xx')
    ax3.set_ylabel('sigma_yy')
    ax3.grid()
    ax3.legend(bbox_to_anchor=(1.5, 1.0), loc='upper right')
    
    canvas.draw()
    canvas.get_tk_widget().grid(row = 0, column=15)
    #canvas.get_tk_widget().grid(row =5, column=0)

if __name__ == '__main__':
             
    root = Tk() 
    root.title('Forward Computation of R values')
    root.geometry("1000x1000")
    #root.geometry("500x500")
    
    count = 0
    
    content = ttk.Frame(root)
    frame1 = ttk.Frame(content, borderwidth=5, relief="ridge", width=700, height=700)
    #frame2 = ttk.Frame(content, borderwidth=5, relief="ridge", width=500, height=450)
    frame3 = ttk.Frame(content, borderwidth=5, relief="ridge", width=700, height=700)
    
    content.grid(column=0, row=0)
    
    frame1.grid(column=0, row=0, columnspan=14, rowspan=14)
    #frame2.grid(column=0, row=15, columnspan=14, rowspan=14)
    frame3.grid(column=15, row=0, columnspan=30, rowspan=30)
    
    fig = plt.figure()
    fig.set_figheight(15)
    fig.set_figwidth(15)
         
    ax1 = plt.subplot2grid(shape=(10, 10), loc=(0, 0), rowspan = 3, colspan=8)
    ax3 = plt.subplot2grid(shape=(10, 10), loc=(4, 1), rowspan=5, colspan=5)
    
    ax2 = ax1.twinx()
    canvas = FigureCanvasTkAgg(fig, master = frame3)
    
    toolbarFrame = Frame(master=frame1)
    toolbarFrame.grid(row=25,column=0)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    
    lng = widgets(frame1, ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'e', 'select'])
    
    root.mainloop()