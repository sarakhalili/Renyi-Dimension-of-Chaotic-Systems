def D_alpha(signal,alpha,L0,scaling_factor,iteration,DimensionofSystem,plotEntropyRadiusOption):
# % This function calculates Renyi generalized Dimension 
# % Author : Sarakhalili
# %------------------------------------------------------------------------------------                                                                      
# %   Input parameters:                                                     
# %       - signal:                     Input signal (There must be enough data inputs otherwise this function would not work properly)                       
# %       - alpha:                      the order of the entropy     
# %       - L0:                         Length of the initial box
# %       - scaling_factor:             the scaling factor which shrinks the size of boxes at each iteration
# %       - iteration:                  number of times scaling factor applies
# %       - DimensionofSignal:          Defines the dimension of system. should be "1D" , "2D" or "3D"
# %       - plotEntropyRadiusOption:    if true plots the logEntropy_logRadius graph
# %   Output:                                                   
# %       - D_alpha:      Renyi generalized Dimension in alpha             
# % -------------------------------------------------------------------------------

  import math
  import numpy as np
  import matplotlib.pyplot as plt
  from scipy.stats import linregress

#=============================================================================================================================================================                                                         
  N=len(signal);
  S = []
  En = []
#=============================================================================================================================================================
  if DimensionofSystem=="1D":
      for i in range(iteration):                                                                  # this loop calculates the Renyi entropy each radius ^ -scaling factor
        p = []
        radius = L0*math.pow(scaling_factor,-i);                                                  # the length of boxes in each dimension
        
        numberofboxes= ((np.amax(signal)-np.amin(signal))/radius);                                # number of boxes created in each dimension

        if math.floor(numberofboxes)!=numberofboxes:                                              # this if checks if the boxes fit the length of signal perfectly in x dimension 
            numberofboxes=numberofboxes+1                                                         # if not adds one to cover the last portion of the signal
        
        numberofboxes=math.floor(numberofboxes)

        for k in range(numberofboxes):                                                            # calculates the probability of each box according to the signal distribution
            m=len(np.argwhere( np.logical_and(((k*radius)+np.amin(signal))<=signal,signal<(((k+1)*radius)+np.amin(signal))) ))
            p.append(m/N)      
        p=np.array(p)                                                                             # deleting p=0 since it has no impact on Entropy
        p=p.flatten()
        p=p[p != 0]
        # print(p)       
        
        if alpha==1:                                                                              # calculates information dimension seperately because Renyi dimension is not defined for alpha=1
            RenyiEn=-sum(p*np.log2(p))
        else:
            RenyiEn=math.log2(sum(np.power(p,alpha)))/(1-alpha);                                  # Renyi dimension

        En.append(RenyiEn)
        S.append(-math.log2(radius))

#===========================================================================================================================================================
  if DimensionofSystem=="2D": 
      x=signal[:,0]; y=signal[:,1]
      for i in range(iteration):                                                                  # this loop calculates the Renyi entropy each radius ^ -scaling factor
        p = []
        radius = L0*math.pow(scaling_factor,-i);                                                  # the length of boxes in each dimension
        print("Number of iteration is : ",i,end="  ")
        
        numberofboxesx= ((np.amax(x)-np.amin(x))/radius);                                         # number of boxes created in x dimension
        numberofboxesy= ((np.amax(y)-np.amin(y))/radius);                                         # number of boxes created in y dimension

        if math.floor(numberofboxesx)!=numberofboxesx:                                            # this if checks if the boxes fit the length of signal perfectly in x dimension 
            numberofboxesx=numberofboxesx+1                                                       # if not adds one to cover the last portion of the signal

        if math.floor(numberofboxesy)!=numberofboxesy:                                            # this if checks if the boxes fit the length of signal perfectly in y dimension 
            numberofboxesy=numberofboxesy+1                                                       # if not adds one to cover the last portion of the signal    


        numberofboxesx=math.floor(numberofboxesx); numberofboxesy=math.floor(numberofboxesy)

        for kx in range(numberofboxesx):                                                          # calculates the probability of each box according to the signal distribution
            condition1 = ( np.argwhere( np.logical_and( ((kx*radius)+np.amin(x))<=x , x<(((kx+1)*radius)+np.amin(x)) ) ) ).flatten()
            for ky in range(numberofboxesy):        
              condition2 = ( np.argwhere( np.logical_and( ((ky*radius)+np.amin(y))<=y , y<(((ky+1)*radius)+np.amin(y)) ) ) ).flatten()    
              intersection = 0
              for c1 in condition1:
                  for c2 in condition2:
                      if c1 == c2:
                          intersection += 1     
              m=intersection/N
              p.append(m)
        p=np.array(p)                                                                             # deleting p=0 since it has no impact on Entropy
        p=p.flatten()
        p=p[p != 0]
        # print(p)       
        
        if alpha==1:                                                                              # calculates information dimension seperately because Renyi dimension is not defined for alpha=1
            RenyiEn=-sum(p*np.log2(p))
        else:
            RenyiEn=math.log2(sum(np.power(p,alpha)))/(1-alpha);                                  # Renyi dimension

        En.append(RenyiEn)
        S.append(-math.log2(radius))

#===========================================================================================================================================================
  if DimensionofSystem=="3D":
      x=signal[:,0]; y=signal[:,1]; z=signal[:,2] 
      for i in range(iteration):                                                               # this loop calculates the Renyi entropy each radius ^ -scaling factor
        p = []
        radius = L0*math.pow(scaling_factor,-i);                                               # the length of boxes in each dimension
        
        # print("Number of iteration is : ",i,end="  ")
        
        numberofboxesx= ((np.amax(x)-np.amin(x))/radius);                                      # number of boxes created in x dimension
        numberofboxesy= ((np.amax(y)-np.amin(y))/radius);                                      # number of boxes created in y dimension
        numberofboxesz= ((np.amax(z)-np.amin(z))/radius);                                      # number of boxes created in y dimension

        if math.floor(numberofboxesx)!=numberofboxesx:                                         # this if checks if the boxes fit the length of signal perfectly in x dimension 
            numberofboxesx=numberofboxesx+1                                                    # if not adds one to cover the last portion of the signal

        if math.floor(numberofboxesy)!=numberofboxesy:                                         # this if checks if the boxes fit the length of signal perfectly in y dimension 
            numberofboxesy=numberofboxesy+1                                                    # if not adds one to cover the last portion of the signal    

        if math.floor(numberofboxesz)!=numberofboxesz:                                         # this if checks if the boxes fit the length of signal perfectly in z dimension 
            numberofboxesz=numberofboxesz+1                                                    # if not adds one to cover the last portion of the signal    


        numberofboxesx=math.floor(numberofboxesx); numberofboxesy=math.floor(numberofboxesy); numberofboxesz=math.floor(numberofboxesz)

        for kx in range(numberofboxesx):                                                        # calculates the probability of each box according to the signal distribution
            condition1 = ( np.argwhere( np.logical_and( ((kx*radius)+np.amin(x))<=x , x<(((kx+1)*radius)+np.amin(x)) ) ) ).flatten()
            for ky in range(numberofboxesy): 
                condition2 = ( np.argwhere( np.logical_and( ((ky*radius)+np.amin(y))<=y , y<(((ky+1)*radius)+np.amin(y)) ) ) ).flatten()

                intersection = 0; cc=[]
                for c1 in condition1:
                    if c1 in condition2:
                        cc.append(c1)
                for kz in range(numberofboxesz):                     
                    condition3 = ( np.argwhere( np.logical_and( ((kz*radius)+np.amin(z))<=z , z<(((kz+1)*radius)+np.amin(z)) ) ) ).flatten()     
                     
                    for cc1 in cc:
                        if cc1 in condition3:
                            intersection += 1     
                    m=intersection/N
                    p.append(m)

        p=np.array(p)                                                                             # deleting p=0 since it has no impact on Entropy
        p=p.flatten()
        p=p[p != 0]    
        
        if alpha==1:                                                                              # calculates information dimension seperately because Renyi dimension is not defined for alpha=1
            RenyiEn=-sum(p*np.log2(p))
        else:
            RenyiEn=math.log2(sum(np.power(p,alpha)))/(1-alpha);                                  # Renyi dimension

        En.append(RenyiEn)
        S.append(-math.log2(radius))
#===========================================================================================================================================================

  slope, intercept, r_value, p_value, std_err = linregress(S, En)                                 # fits a line to the entropy_radius log_log graph                                              
  D_alpha=slope 
  if plotEntropyRadiusOption==True:                                                               # optional: plots the logEntropy_logRadius graph
    xtest=np.linspace(np.amin(S),np.amax(S))
    ytest=slope*xtest+intercept

    plt.plot(S,En,'ro', S, En, 'k',xtest,ytest,'b')        
    plt.title('Entrophy - Radius \n q = {} , Dq = {}'.format(alpha,D_alpha))
    plt.ylabel('Log(Renyi Entrophy)')
    plt.xlabel('Log(Radius)')
    plt.grid(True)
    # plt.savefig('Entrophy - Radius , alpha = {}.png'.format(alpha), dpi=200)
    plt.show() 
  return D_alpha                                                                                  # the slope of the fitted line is the dimension                    