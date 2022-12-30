from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
import sys
import math
#import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from functions_v6 import get_sigVecInRollCS_from_sigPhi, yld2000_derivative, yld2000_EqvStr # , get_R_vs_Theta_Rb_values, get_UniaxYieldStress_vs_Theta_values, get_stress_points_in_Yieldlocus

class widgets(Frame):
    
    def __init__(self, window, picks=[]):
       
      Frame.__init__(self, window)
      
      self.inputs = []
      self.para = ['a2','a3','a4','a5','a6','a7','a8', 'e']
      lf = LabelFrame(window, text='Inputs')
      lf.grid(column=0, row=0)
      value_picks = [[1.0, 1.2, 3], [6,8,3]]
      
      label = Label(lf, text='a1', width=6)
      label.grid(row=0, column=0)      
      var = IntVar()
      var.set(1)
      self.a1 = Entry(lf, textvariable=var, width=6)
      self.a1.grid(row=0, column=1)
      """
      label = Label(lf, text='e_min', width=6)
      label.grid(row=0, column=2)      
      var = IntVar()
      var.set(6)
      self.e_min = Entry(lf, textvariable=var, width=6)
      self.e_min.grid(row=0, column=3)
      
      label = Label(lf, text='e_max', width=6)
      label.grid(row=0, column=4)     
      var = IntVar()
      var.set(8)
      self.e_max = Entry(lf, textvariable=var, width=6)
      self.e_max.grid(row=0, column=5)
      """     
      for j, para in enumerate(self.para): 
          label = Label(lf, text=self.para[j], width=8)
          label.grid(row=0, column=j+3)
                    
          alpha_range  = []
          for i in range(len(picks)): 
              
              label = Label(lf, text=picks[i], width=6)
              label.grid(row=i+1, column=2)
                                
              var = IntVar()
              if para== 'e':
                  var.set(value_picks[1][i])
              else:                  
                  var.set(value_picks[0][i])
                           
              ent = Entry(lf, textvariable=var, width=8)
              ent.grid(row=i+1, column=j+3) 
              
              alpha_range.append(ent)
              
          self.inputs.append(alpha_range)
              
          Button(lf, text='Get the plot!', command=allstates).grid(row=12,column=1)
           
    def state(self): 
              
      a1 = float(self.a1.get())
      #e_min = float(self.e_min.get())
      #e_max = float(self.e_max.get())
      #exponent_values = [float(self.e_min.get()), float(self.e_max.get())]
      
      input_data = []
      for i, alpha_input in enumerate(self.inputs):
          #alpha_values_list = []
          alpha_input = list(map((lambda ent: float(ent.get())), alpha_input))
          #print(alpha_input)
          n_step = alpha_input[-1]
          alpha_values = alpha_input[:-1]
          
          step_alpha = (alpha_values[1]-alpha_values[0]) / (n_step - 1)
          #print(step_alpha)
          if n_step > 2:              
              for j in range(int(n_step-2)): # 2 accounting for the first and the last points                     
                alpha_values = np.insert(alpha_values, j+1, alpha_values[j]+step_alpha, axis=0)
                
          input_data.append(alpha_values)
              
      return input_data, a1, n_step
                   
def allstates():
     
    global count, canvas, fig, ax1, ax2, ax3, toolbar
    
    count +=1 
    
    if count != 1:
        ax1.cla()
        ax2.cla()
        ax3.cla()
        ax4.cla()
        ax5.cla()
        ax6.cla()
        ax7.cla()
        ax8.cla()
    
    input_data, a1, n_step = lng.state()  
       
    print(np.shape(input_data))
    print(input_data)
    
    # R_0_values = []  # or use [0] * inputs[4]**8   OR use numpy arrays (faster and smaller in size)
    # R_45_values = []
    # R_90_values = []
    # R_biax_values = []
    # S_0_values = []
    # S_45_values = []
    # S_90_values = []
    # S_biax_values = []
    
    no_of_calcns = int(n_step)**8
    R_0_values = np.zeros(no_of_calcns).tolist()  # or use [0] * inputs[4]**8   OR use numpy arrays (faster and smaller in size)
    R_45_values = np.zeros(no_of_calcns).tolist()
    R_90_values = np.zeros(no_of_calcns).tolist()
    R_biax_values = np.zeros(no_of_calcns).tolist()
    S_0_values = np.zeros(no_of_calcns).tolist()
    S_45_values = np.zeros(no_of_calcns).tolist()
    S_90_values = np.zeros(no_of_calcns).tolist()
    S_biax_values = np.zeros(no_of_calcns).tolist()
    
    no_of_calcns = 0
    
    
    for a in input_data[7]: # this is the exponent of the yield locus
        for a2 in input_data[0]: 
            for a3 in input_data[1]:
                for a4 in input_data[2]:
                    for a5 in input_data[3]:
                        for a6 in input_data[4]:
                            for a7 in input_data[5]:
                                for a8 in input_data[6]:
                                    
                                    alphaOrj = [a1,a2,a3,a4,a5,a6,a7,a8] 
                                    factor = ( (abs((2.0*alphaOrj[0]+alphaOrj[1])/3.0))**a + (abs((2.0*alphaOrj[2]-2.0*alphaOrj[3])/3.0))**a + (abs((4.0*alphaOrj[4]-alphaOrj[5])/3.0))**a  ) / 2.0
                                    factor = factor ** (1.0/a)
                                    alpha = [i/factor for i in alphaOrj]  # here the normalized values
                                    
                                    YieldStress = 1
                                    angleDeg = 0.0
                                    sig_vec = get_sigVecInRollCS_from_sigPhi(YieldStress, angleDeg) # yield stress and angle in rolling direction
                                    dPhi_dsig = yld2000_derivative(alpha,a,sig_vec)
                                    rad = angleDeg/180*math.pi
                                    R_0 = -(math.sin(rad)**2*dPhi_dsig[0][0] - math.sin(2*rad)*dPhi_dsig[0][2] + math.cos(rad)**2*dPhi_dsig[0][1]) / (dPhi_dsig[0][0] + dPhi_dsig[0][1])
                                    S_0 = YieldStress / yld2000_EqvStr(alpha,a,get_sigVecInRollCS_from_sigPhi(YieldStress, angleDeg))
                                    
                                    # if (-10 > R_0 or R_0 > 10):
                                    #     print("R0 = ", R_0, "a,alpha=", a, alphaOrj)


                                    angleDeg = 45.0
                                    sig_vec = get_sigVecInRollCS_from_sigPhi(YieldStress, angleDeg) # yield stress and angle in rolling direction
                                    dPhi_dsig = yld2000_derivative(alpha,a,sig_vec)
                                    rad = angleDeg/180*math.pi
                                    R_45 = -(math.sin(rad)**2*dPhi_dsig[0][0] - math.sin(2*rad)*dPhi_dsig[0][2] + math.cos(rad)**2*dPhi_dsig[0][1]) / (dPhi_dsig[0][0] + dPhi_dsig[0][1])            
                                    S_45 = YieldStress / yld2000_EqvStr(alpha,a,get_sigVecInRollCS_from_sigPhi(YieldStress, angleDeg))
 
                                    # if (-10 > R_45 or R_45 > 10):
                                    #     print("R45 = ", R_45, "a,alpha=", a, alphaOrj)


                                    angleDeg = 90.0
                                    sig_vec = get_sigVecInRollCS_from_sigPhi(YieldStress, angleDeg) # yield stress and angle in rolling direction
                                    dPhi_dsig = yld2000_derivative(alpha,a,sig_vec)
                                    rad = angleDeg/180*math.pi
                                    R_90 = -(math.sin(rad)**2*dPhi_dsig[0][0] - math.sin(2*rad)*dPhi_dsig[0][2] + math.cos(rad)**2*dPhi_dsig[0][1]) / (dPhi_dsig[0][0] + dPhi_dsig[0][1])
                                    S_90 = YieldStress / yld2000_EqvStr(alpha,a,get_sigVecInRollCS_from_sigPhi(YieldStress, angleDeg))
  
                                    # if (-10 > R_90 or R_90 > 10):
                                    #     print("R90 = ", R_90, "a,alpha=", a, alphaOrj)

                                    sig_vec = np.array([[1.0],[1.0],[0]])
                                    dPhi_dsig_b = yld2000_derivative(alpha,a,sig_vec)
                                    R_biax = dPhi_dsig_b[0][1] / dPhi_dsig_b[0][0]
                                    S_biax = 0.7071067811865476 / yld2000_EqvStr(alpha,a,[0.7071067811865476, 0.7071067811865476, 0.0])
                                    
                                    if (-10 > R_biax or R_biax > 9000):
                                        print("R_biax = ", R_biax, "a,alpha=", a, alphaOrj)
                                   
                                    
                                    R_0_values[no_of_calcns] = R_0  
                                    R_45_values[no_of_calcns] =R_45
                                    R_90_values[no_of_calcns] =R_90
                                    R_biax_values[no_of_calcns] =R_biax
                                    S_0_values[no_of_calcns] =S_0
                                    S_45_values[no_of_calcns] =S_45
                                    S_90_values[no_of_calcns] =S_90
                                    S_biax_values[no_of_calcns] =S_biax
                                    
                                    no_of_calcns += 1
                                       
                                    # R, R_biax = get_R_vs_Theta_Rb_values(alpha, a)
                                    # S = get_UniaxYieldStress_vs_Theta_values(alpha, a, Y_0=1.0, NormalizeAlpha=False)
                                                                           
                                    # R_0, R_45, R_90 = R[0][0], R[0][45], R[0][90]
                                    # S_0, S_45, S_90, S_biax = S[0][2], S[0][4], S[0][3], S[0][5]
                                                                          
                                    # R_0_values.append(R_0)
                                    # R_45_values.append(R_45)
                                    # R_90_values.append(R_90)
                                    # R_biax_values.append(R_biax)
                                    # S_0_values.append(S_0)
                                    # S_45_values.append(S_45)
                                    # S_90_values.append(S_90)
                                    # S_biax_values.append(S_biax)

    
    # print(type(S_biax_values))
    # print(np.shape(S_biax_values))
    
    
    x = list(range(1,no_of_calcns+1))
    ax1.scatter(x, S_0_values)
    ax2.scatter(x, S_45_values)
    ax3.scatter(x, S_90_values)
    ax4.scatter(x, S_biax_values)
    ax5.scatter(x, R_0_values)
    ax6.scatter(x, R_45_values)
    ax7.scatter(x, R_90_values)
    ax8.scatter(x, R_biax_values)
    
    min_S_0, max_S_0   = min(S_0_values), max(S_0_values)
    min_S_45, max_S_45 = min(S_45_values), max(S_45_values)
    min_S_90, max_S_90 = min(S_90_values), max(S_90_values)
    min_S_biax, max_S_biax = min(S_biax_values), max(S_biax_values)
    min_R_0, max_R_0   = min(R_0_values), max(R_0_values)
    min_R_45, max_R_45 = min(R_45_values), max(R_45_values)
    min_R_90, max_R_90 = min(R_90_values), max(R_90_values)
    min_R_biax, max_R_biax = min(R_biax_values), max(R_biax_values)
    
    ax1.set_ylabel('S_0', fontweight='bold')
    ax1.set_xlabel('no. of calculations', fontweight='bold') 
    ax2.set_ylabel('S_45', fontweight='bold')
    ax2.set_xlabel('no. of calculations', fontweight='bold')  
    ax3.set_ylabel('S_90', fontweight='bold')
    ax3.set_xlabel('no. of calculations', fontweight='bold') 
    ax4.set_ylabel('S_biax', fontweight='bold')
    ax4.set_xlabel('no. of calculations', fontweight='bold') 
    ax5.set_ylabel('R_0', fontweight='bold')
    ax5.set_xlabel('no. of calculations', fontweight='bold')
    ax6.set_ylabel('R_45', fontweight='bold')
    ax6.set_xlabel('no. of calculations', fontweight='bold') 
    ax7.set_ylabel('R_90', fontweight='bold')
    ax7.set_xlabel('no. of calculations', fontweight='bold') 
    ax8.set_ylabel('R_biax', fontweight='bold')
    ax8.set_xlabel('no. of calculations', fontweight='bold') 
    
    ax1.axhline(y=max_S_0, color='red', label=f'max={max_S_0}')
    ax1.axhline(y=min_S_0, color='green', label=f'min={min_S_0}')
    ax1.legend(loc='upper right')
    ax2.axhline(y=max_S_45, color='red', label=f'max={max_S_45}')
    ax2.axhline(y=min_S_45, color='green', label=f'min={min_S_45}')
    ax2.legend(loc='upper right')
    ax3.axhline(y=max_S_90, color='red', label=f'max={max_S_90}')
    ax3.axhline(y=min_S_90, color='green', label=f'min={min_S_90}')
    ax3.legend(loc='upper right')
    ax4.axhline(y=max_S_biax, color='red', label=f'max={max_S_biax}')
    ax4.axhline(y=min_S_biax, color='green', label=f'min={min_S_biax}')
    ax4.legend(loc='upper right')
    ax5.axhline(y=max_R_0, color='red', label=f'max={max_R_0[0]}')
    ax5.axhline(y=min_R_0, color='green', label=f'min={min_R_0[0]}')
    ax5.legend(loc='upper right')
    ax6.axhline(y=max_R_45, color='red', label=f'max={max_R_45[0]}')
    ax6.axhline(y=min_R_45, color='green', label=f'min={min_R_45[0]}')
    ax6.legend(loc='upper right')
    ax7.axhline(y=max_R_90, color='red', label=f'max={max_R_90[0]}')
    ax7.axhline(y=min_R_90, color='green', label=f'min={min_R_90[0]}')
    ax7.legend(loc='upper right')
    ax8.axhline(y=max_R_biax, color='red', label=f'max={max_R_biax[0]}')
    ax8.axhline(y=min_R_biax, color='green', label=f'min={min_R_biax[0]}')
    ax8.legend(loc='upper right')
        
    #ax1.tick_params(axis='both')
    #ax2.tick_params(axis='both')
    #ax3.tick_params(axis='both')

    canvas.draw()
    canvas.get_tk_widget().grid(row = 10, column=0)
    
    ##################SAVE INPUT AND STRESS AND R VALUES#################################
    """
    filename =  f'stress_R_values_a_min={inputs[0]}_amax={inputs[1]}_alpha_min={inputs[2]}_alpha_max={inputs[3]}_NPoints={inputs[4]}.h5'
    F = h5py.File(filename, 'a')
    	
    group1 = F.require_group('input_data')    
    dset_11 = group1.require_dataset('input_data', data= input_data, shape=np.shape(input_data), dtype=np.float64, compression='gzip')
    	
    group2 = F.require_group('output_data')    
    dset_21 = group2.require_dataset('R_0_values', data= R_0_values, shape=np.shape(R_0_values), dtype=np.float64, compression='gzip')
    dset_22 = group2.require_dataset('R_45_values', data= R_45_values, shape=np.shape(R_45_values), dtype=np.float64, compression='gzip')
    dset_23 = group2.require_dataset('R_90_values', data= R_90_values, shape=np.shape(R_90_values), dtype=np.float64, compression='gzip')
    dset_24 = group2.require_dataset('R_biax_values', data= R_biax_values, shape=np.shape(R_biax_values), dtype=np.float64, compression='gzip')
    dset_25 = group2.require_dataset('S_0_values', data= S_0_values, shape=np.shape(S_0_values), dtype=np.float64, compression='gzip')
    dset_26 = group2.require_dataset('S_45_values', data= S_45_values, shape=np.shape(S_45_values), dtype=np.float64, compression='gzip')
    dset_27 = group2.require_dataset('S_90_values', data= S_90_values, shape=np.shape(S_90_values), dtype=np.float64, compression='gzip')
    dset_28 = group2.require_dataset('S_biax_values', data= S_biax_values, shape=np.shape(S_biax_values), dtype=np.float64, compression='gzip')    
    """								
    ###########################################################################
    
    R_0_values.clear()
    R_45_values.clear()
    R_90_values.clear()
    R_biax_values.clear()
    S_0_values.clear()
    S_45_values.clear()
    S_90_values.clear()
    S_biax_values.clear()
    
    print('plotted!')
                                
if __name__ == '__main__':
             
    root = Tk() 
    root.title('Forward Computation')
    root.geometry("1280x800")
    #root.geometry("500x500")
    
    count = 0
    
    content = ttk.Frame(root)
    frame1 = ttk.Frame(content, borderwidth=5, relief="ridge")
    frame3 = ttk.Frame(content, borderwidth=5, relief="ridge")
    #width and height here sets the dimension of the frame , width=1250, height=550
    content.grid(column=0, row=0)
    
    frame1.grid(column=0, row=0)
    frame3.grid(column=0, row=10, columnspan=70, rowspan=70)
    
    fig = plt.figure()
    fig.set_figheight(8) #sets the height and width of the white space inside frame3
    fig.set_figwidth(18)
         
    ax1 = plt.subplot2grid(shape=(2,4), loc=(0, 0), rowspan = 1, colspan = 1) #shape=Shape of grid in which to place axis
    ax2 = plt.subplot2grid(shape=(2,4), loc=(0, 1), rowspan = 1, colspan = 1)
    ax3 = plt.subplot2grid(shape=(2,4), loc=(0, 2), rowspan = 1, colspan = 1)
    ax4 = plt.subplot2grid(shape=(2,4), loc=(0, 3), rowspan = 1, colspan = 1)
    ax5 = plt.subplot2grid(shape=(2,4), loc=(1, 0), rowspan = 1, colspan = 1)
    ax6 = plt.subplot2grid(shape=(2,4), loc=(1, 1), rowspan = 1, colspan = 1)
    ax7 = plt.subplot2grid(shape=(2,4), loc=(1, 2), rowspan = 1, colspan = 1)
    ax8 = plt.subplot2grid(shape=(2,4), loc=(1, 3), rowspan = 1, colspan = 1)
    
    #fig.tight_layout()
    
    canvas = FigureCanvasTkAgg(fig, master = frame3)
    
    toolbarFrame = Frame(master=frame1)
    toolbarFrame.grid(row=25,column=0)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    
    lng = widgets(frame1, ['min', 'max', 'NPoints'])
    
    root.mainloop()