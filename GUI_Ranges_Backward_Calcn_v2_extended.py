from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
import sys
import math
import time
import os
import numpy as np
import matplotlib.pyplot as plt
from functions_v6 import calculate_optimized_alphas, get_stress_points_in_Yieldlocus, get_R_vs_Theta_Rb_values

class widgets(Frame):
    
    def __init__(self, window, picks=[]):
       
      Frame.__init__(self, window)
      
      self.inputs = []
      self.stress_and_R = ['S45','S90','Sbiax','R00','R45','R90','Rbiax', 'e']
      lf = LabelFrame(window, text='Inputs')
      lf.grid(column=0, row=0)
      value_picks = [[0.85, 1.15, 2], [6,8,2]]
      
      label = Label(lf, text='S00', width=6)
      label.grid(row=0, column=0)      
      var = IntVar()
      var.set(1)
      self.S00 = Entry(lf, textvariable=var, width=6)
      self.S00.grid(row=0, column=1)
      
      for j, para in enumerate(self.stress_and_R): 
          label = Label(lf, text=self.stress_and_R[j], width=8)
          label.grid(row=0, column=j+3)
                    
          stress_and_R_range  = []
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
              
              stress_and_R_range.append(ent)
              
          self.inputs.append(stress_and_R_range)
              
          Button(lf, text='Get the plot!', command=allstates).grid(row=12,column=1)
           
    def state(self): 
              
      S00 = float(self.S00.get())
      #e_min = float(self.e_min.get())
      #e_max = float(self.e_max.get())
      #exponent_values = [float(self.e_min.get()), float(self.e_max.get())]
      
      input_data = []
      for i, stress_and_R_input in enumerate(self.inputs):
          #alpha_values_list = []
          stress_and_R_input = list(map((lambda ent: float(ent.get())), stress_and_R_input))
          #print(alpha_input)
          n_step = stress_and_R_input[-1]
          stress_and_R_values = stress_and_R_input[:-1]
          
          step = (stress_and_R_values[1]-stress_and_R_values[0]) / (n_step - 1)
          #print(step_alpha)
          if n_step > 2:              
              for j in range(int(n_step-2)): # 2 accounting for the first and the last points                     
                stress_and_R_values = np.insert(stress_and_R_values, j+1, stress_and_R_values[j]+step, axis=0)
                
          input_data.append(stress_and_R_values)
              
      return input_data, S00, n_step
                   
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
    
    input_data, S00, n_step = lng.state()  
       
    print(np.shape(input_data))
    print(input_data)
    
    no_of_calcns = int(n_step)**8
    alpha1_values = np.zeros(no_of_calcns).tolist()  # or use [0] * inputs[4]**8   OR use numpy arrays (faster and smaller in size)
    alpha2_values = np.zeros(no_of_calcns).tolist()
    alpha3_values = np.zeros(no_of_calcns).tolist()
    alpha4_values = np.zeros(no_of_calcns).tolist()
    alpha5_values = np.zeros(no_of_calcns).tolist()
    alpha6_values = np.zeros(no_of_calcns).tolist()
    alpha7_values = np.zeros(no_of_calcns).tolist()
    alpha8_values = np.zeros(no_of_calcns).tolist()
    
    no_of_calcns = 0
    track_calcn_no = []   
    for a in input_data[7]: # this is the exponent of the yield locus
        for S45 in input_data[0]: 
            for S90 in input_data[1]:
                for Sbiax in input_data[2]:
                    for R00 in input_data[3]:
                        for R45 in input_data[4]:
                            for R90 in input_data[5]:
                                for Rbiax in input_data[6]:
                                    
                                    alpha = calculate_optimized_alphas(a,S00,S45,S90,Sbiax,R00,R45,R90,Rbiax)
                                    
                                    S = get_stress_points_in_Yieldlocus(alpha,a,shear=0)  
                                    R,Rb = get_R_vs_Theta_Rb_values(alpha,a)
                                    
                                    given_stress = np.array([[S00, S45, S90, Sbiax]])
                                    given_R = np.array([[R00, R45, R90, Rbiax]])
                                        
                                    calculated_stress = np.array([[S[0][2], S[0][4], S[0][3], S[0][5]]])
                                    calculated_R = np.array([[R[0][0], R[0][45], R[0][90], Rb[0]]])

                                    Absolute_Error_for_stress = (abs(given_stress - calculated_stress))
                                    Absolute_Error_for_R_value = (abs(given_R - calculated_R))
                                    # the lines below perform element by element division, / is a shorthand for numpy.divide
                                    Percentage_Error_for_stress = (abs(given_stress - calculated_stress) / given_stress)*100.0
                                    Percentage_Error_for_R_value = (abs(given_R - calculated_R) / given_R)*100.0
                                                                            
                                    alpha1_values[no_of_calcns] = alpha[0] 
                                    alpha2_values[no_of_calcns] = alpha[1]
                                    alpha3_values[no_of_calcns] = alpha[2]
                                    alpha4_values[no_of_calcns] = alpha[3]
                                    alpha5_values[no_of_calcns] = alpha[4]
                                    alpha6_values[no_of_calcns] = alpha[5]
                                    alpha7_values[no_of_calcns] = alpha[6]
                                    alpha8_values[no_of_calcns] = alpha[7]
                                    
                                    for (pe_stress, pe_R) in zip(Percentage_Error_for_stress[0], Percentage_Error_for_R_value[0]):
                                        if (pe_stress > 1) or (pe_R > 1): 
                                            track_calcn_no.append(no_of_calcns)
                                            break
                                            # if no_of_calcns not in track_calcn_no:
                                            #     track_calcn_no.append(no_of_calcns)
                                                
                                    no_of_calcns += 1
                                    
                                    #break
                                       
    x = list(range(no_of_calcns))
    print("total number of simulations:", no_of_calcns)
    
    min_alpha1, max_alpha1 = min(alpha1_values), max(alpha1_values)
    min_alpha2, max_alpha2 = min(alpha2_values), max(alpha2_values)
    min_alpha3, max_alpha3 = min(alpha3_values), max(alpha3_values)
    min_alpha4, max_alpha4 = min(alpha4_values), max(alpha4_values)
    min_alpha5, max_alpha5 = min(alpha5_values), max(alpha5_values)
    min_alpha6, max_alpha6 = min(alpha6_values), max(alpha6_values)
    min_alpha7, max_alpha7 = min(alpha7_values), max(alpha7_values)
    min_alpha8, max_alpha8 = min(alpha8_values), max(alpha8_values)
    
    print("red points:", track_calcn_no)
    
    for idx in x:
        if idx in track_calcn_no:
            ax1.scatter(idx, alpha1_values[idx], color='red', s=10)
            ax2.scatter(idx, alpha2_values[idx], color='red', s=10)
            ax3.scatter(idx, alpha3_values[idx], color='red', s=10)
            ax4.scatter(idx, alpha4_values[idx], color='red', s=10)
            ax5.scatter(idx, alpha5_values[idx], color='red', s=10)
            ax6.scatter(idx, alpha6_values[idx], color='red', s=10)
            ax7.scatter(idx, alpha7_values[idx], color='red', s=10)
            ax8.scatter(idx, alpha8_values[idx], color='red', s=10)
            
        else:
            ax1.scatter(idx, alpha1_values[idx], color='green', s=10)
            ax2.scatter(idx, alpha2_values[idx], color='green', s=10)
            ax3.scatter(idx, alpha3_values[idx], color='green', s=10)
            ax4.scatter(idx, alpha4_values[idx], color='green', s=10)
            ax5.scatter(idx, alpha5_values[idx], color='green', s=10)
            ax6.scatter(idx, alpha6_values[idx], color='green', s=10)
            ax7.scatter(idx, alpha7_values[idx], color='green', s=10)
            ax8.scatter(idx, alpha8_values[idx], color='green', s=10)
        
    ax1.set_ylabel('alpha1', fontweight='bold')
    ax1.set_xlabel('no. of calculations', fontweight='bold') 
    ax2.set_ylabel('alpha2', fontweight='bold')
    ax2.set_xlabel('no. of calculations', fontweight='bold')  
    ax3.set_ylabel('alpha3', fontweight='bold')
    ax3.set_xlabel('no. of calculations', fontweight='bold') 
    ax4.set_ylabel('alpha4', fontweight='bold')
    ax4.set_xlabel('no. of calculations', fontweight='bold') 
    ax5.set_ylabel('alpha5', fontweight='bold')
    ax5.set_xlabel('no. of calculations', fontweight='bold')
    ax6.set_ylabel('alpha6', fontweight='bold')
    ax6.set_xlabel('no. of calculations', fontweight='bold') 
    ax7.set_ylabel('alpha7', fontweight='bold')
    ax7.set_xlabel('no. of calculations', fontweight='bold') 
    ax8.set_ylabel('alpha8', fontweight='bold')
    ax8.set_xlabel('no. of calculations', fontweight='bold') 
    
    ax1.axhline(y=max_alpha1, color='red', label=f'max={max_alpha1}')
    ax1.axhline(y=min_alpha1, color='green', label=f'min={min_alpha1}')
    ax1.legend(loc='upper right')
    ax2.axhline(y=max_alpha2, color='red', label=f'max={max_alpha2}')
    ax2.axhline(y=min_alpha2, color='green', label=f'min={min_alpha2}')
    ax2.legend(loc='upper right')
    ax3.axhline(y=max_alpha3, color='red', label=f'max={max_alpha3}')
    ax3.axhline(y=min_alpha3, color='green', label=f'min={min_alpha3}')
    ax3.legend(loc='upper right')
    ax4.axhline(y=max_alpha4, color='red', label=f'max={max_alpha4}')
    ax4.axhline(y=min_alpha4, color='green', label=f'min={min_alpha4}')
    ax4.legend(loc='upper right')
    ax5.axhline(y=max_alpha5, color='red', label=f'max={max_alpha5}')
    ax5.axhline(y=min_alpha5, color='green', label=f'min={min_alpha5}')
    ax5.legend(loc='upper right')
    ax6.axhline(y=max_alpha6, color='red', label=f'max={max_alpha6}')
    ax6.axhline(y=min_alpha6, color='green', label=f'min={min_alpha6}')
    ax6.legend(loc='upper right')
    ax7.axhline(y=max_alpha7, color='red', label=f'max={max_alpha7}')
    ax7.axhline(y=min_alpha7, color='green', label=f'min={min_alpha7}')
    ax7.legend(loc='upper right')
    ax8.axhline(y=max_alpha8, color='red', label=f'max={max_alpha8}')
    ax8.axhline(y=min_alpha8, color='green', label=f'min={min_alpha8}')
    ax8.legend(loc='upper right')
        
    #ax1.tick_params(axis='both')
    #ax2.tick_params(axis='both')
    #ax3.tick_params(axis='both')

    canvas.draw()
    canvas.get_tk_widget().grid(row = 10, column=0)
				
    ###########################################################################
    
    alpha1_values.clear()
    alpha2_values.clear()
    alpha3_values.clear()
    alpha4_values.clear()
    alpha5_values.clear()
    alpha6_values.clear()
    alpha7_values.clear()
    alpha8_values.clear()
    
    print('plotted!')
    
    end_time = time.time()
    total_time = (end_time - start_time)/60
    print(f'Total Time: {(total_time):.4f} mins')
                                  
if __name__ == '__main__':
             
    root = Tk() 
    root.title('Reverse  Computation')
    root.geometry("1280x800")
    
    start_time = time.time()
    
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