from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from PIL import ImageTk, Image
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
import math
import numpy as np
import matplotlib.pyplot as plt
from functions_v6 import get_R_vs_Theta_Rb_values, get_UniaxYieldStress_vs_Theta_values, get_stress_points_in_Yieldlocus

class widgets(Frame):
    
    def __init__(self, window, picks=[]):
       
      Frame.__init__(self, window)
      
      self.curves = []
      self.base_curve = []
      self.curvenames = ['B.C','C1','C2','C3','C4','C5']
      lf = LabelFrame(window, text='Parameters')
      lf.grid(column=0, row=0)
      
      ############ LABEL AND ENTRY FOR DELTA FOR ALPHAS AND EXPONENT ##########
      label = Label(lf, text='delta', width=6)
      label.grid(row=0, column=1)
      
      var = IntVar()
      var.set(0.1)
      self.delta = Entry(lf, textvariable=var, width=6)
      self.delta.grid(row=0, column=2)
      
      """
      label = Label(lf, text='dExp', width=6)
      label.grid(row=0, column=3)
           
      var = IntVar()
      var.set(0)
      self.dExp = Entry(lf, textvariable=var, width=6)
      self.dExp.grid(row=0, column=4) 
      """
      #########################################################################
      
      for j in range(len(self.curvenames)): 
          
          label = Label(lf, text=self.curvenames[j], width=6)
          label.grid(row=j+2, column=0)
          
          self.alpha_checks = []
          for i, pick in enumerate(picks):
              
              if j == 0:
                  label = Label(lf, text=pick, width=6)
                  label.grid(row=j+1, column=i+1)
                  
              if pick != 'select':
                                              
                  var = IntVar()                            
                  if pick == 'e':
                      var.set(8)
                  else:
                      var.set(1)
                      
                  ent = Entry(lf, textvariable=var, width=6)                     
                  ent.grid(row=2, column=i+1)   
                  
                  if j == len(self.curvenames)-1:  ###
                      self.base_curve.append(ent) 
                      
                  #chk = IntVar()
                  #cb = Checkbutton(lf, variable=chk)
              
              if j != len(self.curvenames)-1:
                  
                  chk = IntVar()
                  cb = Checkbutton(lf, variable=chk)
                  
                  cb.grid(row=j+3,column=i+1)
              self.alpha_checks.append(chk)
                       
          self.curves.append(self.alpha_checks)
          
      Button(lf, text='Get the plot!', command=allstates).grid(row=(len(self.curvenames)+3),column=1)  
      Button(lf, text='Normalize the plot!', command=normalize).grid(row=(len(self.curvenames)+3),column=3)
      
    def state(self): 
      
      delta = float(self.delta.get())

      result_for_cvs = []
      result_back_cvs = []
      result_names =  []
      
      for cu, curve in enumerate(self.curves[:-1]): ###

          curve = list(map((lambda chk: chk.get()), curve))
          #print(curve)
          bc_for = list(map((lambda ent: float(ent.get())), self.base_curve))
          bc_back = list(map((lambda ent: float(ent.get())), self.base_curve))
          if curve[-1] ==1:
              for ch, check in enumerate(curve[:-1]):   
                  #print(check)
                  if check ==1:
                      bc_for[ch] = bc_for[ch] + delta
                      bc_back[ch] = bc_back[ch] - delta
                      
              result_for_cvs.append(bc_for)
              result_back_cvs.append(bc_back)
              result_names.append(self.curvenames[cu+1])

      return delta, result_for_cvs, result_back_cvs, result_names

class display(Frame):
    
    def __init__(self, window, values, curve_names, picks=[]):
       
      Frame.__init__(self, window)
      
      lf = LabelFrame(window, text='Values')
      lf.grid(column=0, row=15)
      #lf.grid_propagate(True)
      
      for j, value in enumerate(values):
          
          label = Label(lf, text=curve_names[j], width=10)
          label.grid(row=j+1, column=0)
          
          #self.vars = []
          
          for i, pick in enumerate(picks):
              
              if j == 0:
                  label = Label(lf, text=pick, width=10)
                  label.grid(row=j, column=i+1)
              
              var =  IntVar(value=value[i])
              ent = Entry(lf, textvariable=var, width=10)
              ent.grid(row=j+1, column=i+1)
                                             
def allstates():
     
    global count, canvas, fig, ax1, ax2, ax3, toolbar
    
    count +=1 
    
    if count != 1:
        ax1.cla()
        ax2.cla()
        ax3.cla()
        
    delta, for_curves, back_curves, curve_names = lng.state()  
    print(delta)
    print(for_curves)
    print(back_curves)
    print(curve_names)

    deg = np.array([list(range(0,91))])
    #deg_yl = np.array([list(range(0,361))])
    deg_yl = np.array([list(range(-45,136))])
    print(deg_yl.shape)
    
    values = []    
    for l, (for_curve, back_curve) in enumerate(zip(for_curves, back_curves)):
                
        for_R, for_Rb = get_R_vs_Theta_Rb_values(for_curve[:-1], for_curve[-1])
        for_stress = get_UniaxYieldStress_vs_Theta_values(for_curve[:-1], for_curve[-1], Y_0=1.0, NormalizeAlpha=False)
        for_yl = get_stress_points_in_Yieldlocus(for_curve[:-1], for_curve[-1], shear=0)
                     
        back_R, back_Rb = get_R_vs_Theta_Rb_values(back_curve[:-1], back_curve[-1])
        back_stress = get_UniaxYieldStress_vs_Theta_values(back_curve[:-1], back_curve[-1], Y_0=1.0, NormalizeAlpha=False)
        back_yl = get_stress_points_in_Yieldlocus(back_curve[:-1], back_curve[-1],shear=0)
        
        sens_stress = np.absolute(for_stress - back_stress)/ (2*delta)
        sens_R = np.absolute (for_R - back_R)/ (2*delta)       
        sens_yl = np.sqrt( ((for_yl[0][0][0]-back_yl[0][0][0])**2) +  ((for_yl[0][1][0]-back_yl[0][1][0])**2) )  / (2*delta)         
        sens_yl_315to135 = []
        sens_yl_315to135.extend(sens_yl[315:360])
        sens_yl_315to135.extend(sens_yl[0:136])
        #print(np.shape(sens_yl_315to135))
        
        angle_plane_strain_1_for = for_yl[0][12]
        angle_plane_strain_2_for = for_yl[0][13]
        angle_sigma_biax_for = for_yl[0][14] 
        angle_plane_strain_1_back = back_yl[0][12]
        angle_plane_strain_2_back = back_yl[0][13]
        angle_sigma_biax_back = back_yl[0][14]
        
        ax1.scatter(deg, sens_R, marker='d')      
        ax2.scatter(deg, sens_stress, marker='.', linewidths=0.0)
        ax3.plot(deg_yl[0], sens_yl_315to135,  label = curve_names[l])
        
        v = [angle_plane_strain_1_for, angle_plane_strain_2_for, angle_sigma_biax_for,
             angle_plane_strain_1_back, angle_plane_strain_2_back, angle_sigma_biax_back]        
        values.append(v)
        
    display(frame1, values, curve_names, ['angle_PS1x_f', 'angle_PS2y_f', 'angle_Sbiax_f',
                                         'angle_PS1x_b', 'angle_PS2y_b', 'angle_Sbiax_b' ])  

    ax1.set_ylabel('R-values Sensitivity\u2666 ',fontsize=12)
    ax1.set_xlabel('angle (in degree) wrt RD',fontsize=12)  
    ax1.set_title('Parameter Sensitivity Plot for Lankford Anistropy Parameters',fontsize=12)      
    ax2.set_ylabel('Stress-values Sensitivity\u2022',fontsize=12)
    ax3.set_title('Parameter Sensitivity Plot for Yield Locus',fontsize=12)
    #ax1.grid()
    ax3.set_xlabel('angle (in degree)',fontsize=12)
    ax3.set_ylabel('Yield Locus Sensitivity',fontsize=12)
    #ax3.grid()
    #ax3.axvline(x=0, color='#FF3339', label='sigma_00', ls='--')
    #ax3.axvline(x=90, color='#F615DE', label='sigma_90', ls='--')
    #ax3.axvline(x=angle_plane_strain_1, color='#070707', label=f'sigma_plane_strain_1={angle_plane_strain_1:.4f}', ls='--')
    #ax3.axvline(x=angle_plane_strain_2, color='#60A380', label=f'sigma_plane_strain_2={angle_plane_strain_2:.4f}', ls='--')
    #ax3.axvline(x=angle_sigma_biax, color='#15F669', label=f'sigma_biax={angle_sigma_biax:.4f}', ls='--')
    #ax3.axvline(x=135, color='#F3F615', label='sigma_shear_1', ls='--')
    #ax3.axvline(x=-45, color='#F6AE15', label='sigma_shear_2', ls='--')
    #ax3.axvline(x=45, color='#1521F6', label='sigma_45', ls='--')
    ax3.legend(bbox_to_anchor=(1.6, 1.0), loc='upper right')
    
    ax1.tick_params(axis='both', labelsize=12)
    ax2.tick_params(axis='both', labelsize=12)
    ax3.tick_params(axis='both', labelsize=12)
        
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
        
    delta, for_curves, back_curves, curve_names = lng.state()
    
    def normalize_alpha(alpha,a):
        factor = ( (abs((2.0*alpha[0]+alpha[1])/3.0))**a + (abs((2.0*alpha[2]-2.0*alpha[3])/3.0))**a + (abs((4.0*alpha[4]-alpha[5])/3.0))**a  ) / 2.0
        factor =  round(float(factor ** (1.0/a)),4)
        norm_alpha = [i/factor for i in alpha]  # here the normalized values
        return norm_alpha, a
    
    deg = np.array([list(range(0,91))])
    deg_yl = np.array([list(range(-45,136))])
      
    values = []  
    for l, (for_curve, back_curve) in enumerate(zip(for_curves, back_curves)):
        
        for_alpha, for_a = normalize_alpha(for_curve[:-1], for_curve[-1])
        for_R, for_Rb = get_R_vs_Theta_Rb_values(for_alpha, for_a)
        for_stress = get_UniaxYieldStress_vs_Theta_values(for_alpha, for_a, Y_0=1.0, NormalizeAlpha=False)
        for_yl = get_stress_points_in_Yieldlocus(for_alpha, for_a, shear=0)
        
        back_alpha, back_a = normalize_alpha(back_curve[:-1], back_curve[-1])           
        back_R, back_Rb = get_R_vs_Theta_Rb_values(back_alpha, back_a)
        back_stress = get_UniaxYieldStress_vs_Theta_values(back_alpha, back_a, Y_0=1.0, NormalizeAlpha=False)
        back_yl = get_stress_points_in_Yieldlocus(back_alpha, back_a, shear=0)
        
        sens_stress = np.absolute(for_stress - back_stress)/ (2*delta)
        sens_R = np.absolute (for_R - back_R)/ (2*delta)       
        sens_yl = np.sqrt( ((for_yl[0][0][0]-back_yl[0][0][0])**2) +  ((for_yl[0][1][0]-back_yl[0][1][0])**2) )  / (2*delta) 
        
        sens_yl_315to135 = []
        sens_yl_315to135.extend(sens_yl[315:360])
        sens_yl_315to135.extend(sens_yl[0:136])
        #print(np.shape(sens_yl_315to135))
        
        angle_plane_strain_1_for = for_yl[0][12]
        angle_plane_strain_2_for = for_yl[0][13]
        angle_sigma_biax_for = for_yl[0][14] 
        angle_plane_strain_1_back = back_yl[0][12]
        angle_plane_strain_2_back = back_yl[0][13]
        angle_sigma_biax_back = back_yl[0][14]
        
        ax1.scatter(deg, sens_R, marker='d')      
        ax2.scatter(deg, sens_stress, marker='.', linewidths=0.0)
        ax3.plot(deg_yl[0], sens_yl_315to135,  label = curve_names[l])
        
        v = [angle_plane_strain_1_for, angle_plane_strain_2_for, angle_sigma_biax_for,
             angle_plane_strain_1_back, angle_plane_strain_2_back, angle_sigma_biax_back]        
        values.append(v)
        
    display(frame1, values, curve_names, ['angle_PS1x_f', 'angle_PS2y_f', 'angle_Sbiax_f',
                                         'angle_PS1x_b', 'angle_PS2y_b', 'angle_Sbiax_b' ])     

    ax1.set_ylabel('R-values Sensitivity\u2666 ',fontsize=12)
    ax1.set_xlabel('angle (in degree) wrt RD',fontsize=12)  
    ax1.set_title('Parameter Sensitivity Plot for Lankford Anistropy Parameters',fontsize=12)      
    ax2.set_ylabel('Stress-values Sensitivity\u2022',fontsize=12)
    ax3.set_title('Parameter Sensitivity Plot for Yield Locus',fontsize=12)
    #ax1.grid()
    ax3.set_xlabel('angle (in degree)',fontsize=12)
    ax3.set_ylabel('Yield Locus Sensitivity',fontsize=12)
    #ax3.grid()
    #ax3.axvline(x=0, color='#FF3339', label='sigma_00', ls='--')
    #ax3.axvline(x=90, color='#F615DE', label='sigma_90', ls='--')
    #ax3.axvline(x=angle_plane_strain_1, color='#070707', label=f'sigma_plane_strain_1={angle_plane_strain_1:.4f}', ls='--')
    #ax3.axvline(x=angle_plane_strain_2, color='#60A380', label=f'sigma_plane_strain_2={angle_plane_strain_2:.4f}', ls='--')
    #ax3.axvline(x=angle_sigma_biax, color='#15F669', label=f'sigma_biax={angle_sigma_biax:.4f}', ls='--')
    #ax3.axvline(x=135, color='#F3F615', label='sigma_shear_1', ls='--')
    #ax3.axvline(x=-45, color='#F6AE15', label='sigma_shear_2', ls='--')
    #ax3.axvline(x=45, color='#1521F6', label='sigma_45', ls='--')
    ax3.legend(bbox_to_anchor=(1.6, 1.0), loc='upper right')
    
    ax1.tick_params(axis='both', labelsize=12)
    ax2.tick_params(axis='both', labelsize=12)
    ax3.tick_params(axis='both', labelsize=12)
        
    canvas.draw()
    canvas.get_tk_widget().grid(row = 0, column=15)
    #canvas.get_tk_widget().grid(row =5, column=0)

if __name__ == '__main__':
             
    root = Tk() 
    root.title('Forward Computation')
    root.geometry("1280x800")
    #root.geometry("500x500")
    
    count = 0
    
    content = ttk.Frame(root)
    frame1 = ttk.Frame(content, borderwidth=5, relief="ridge")
    #frame2 = ttk.Frame(content, borderwidth=5, relief="ridge", width=500, height=450)
    frame3 = ttk.Frame(content, borderwidth=5, relief="ridge")
    
    content.grid(column=0, row=0)
    
    frame1.grid(column=0, row=0, columnspan=14, rowspan=14)
    #frame2.grid(column=0, row=15, columnspan=14, rowspan=14)
    frame3.grid(column=15, row=0, columnspan=30, rowspan=30)
    
    fig = plt.figure()
    fig.set_figheight(10)
    fig.set_figwidth(9)
         
    ax1 = plt.subplot2grid(shape=(10, 10), loc=(0, 0), rowspan = 3, colspan=8)
    ax3 = plt.subplot2grid(shape=(10, 10), loc=(5, 0), rowspan=4, colspan=6)
    
    ax2 = ax1.twinx()
    canvas = FigureCanvasTkAgg(fig, master = frame3)
    
    toolbarFrame = Frame(master=frame1)
    toolbarFrame.grid(row=25,column=0)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    
    img = Image.open("sensitivity.jpg")
    img = img.resize((500,300), Image.ANTIALIAS)
    img = ImageTk.PhotoImage(img)
    label = Label(frame1, image = img).grid(row=30,column=0)
    
    lng = widgets(frame1, ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'e', 'select'])
    
    root.mainloop()