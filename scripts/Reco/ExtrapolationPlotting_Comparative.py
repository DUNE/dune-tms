#These will be useful
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#A brief function to turn the input "num" below into a proper label (inserts a comma)
def int_to_string(integer):
    start_string=str(integer)
    if len(start_string)!=3:
        print("I couldn't get those parameters")
    else:
        display_string=start_string[0]+","+start_string[1:]
    return display_string

#Function for making histogram data
#Inputs: Type - the coordinate/position we investigate e.g. UEndZ, num - the choice of parameters (expressed as a 3-digit integer), range_start - start limit for histograms,range_end - end limit for the histograms, colour - the colour to plot the line in, tag - if this parameter choice is "Optimum" or "Sub-optimum"
#Outputs: The plotted histograms, and the handles for the legend in an array
def DataProcessor(Type,num,range_start,range_end,colour,tag): 
    #Read in a single data file of the Type and num
    MyData=np.loadtxt("/uni-mainz.de/homes/bellbenj/Documents/Extrapolation Stuff/Histogram%sData%d.txt"%(Type,num), float)
    #Prepare the data for plotting by treating the read-in file and make a histogram out of it
    if Type[-1]=="Z":
        MyData=MyData[0:len(MyData)-1:2]/10
        plt.hist(MyData[0:-1],bins=60,fill=False,histtype='step',range=(range_start,range_end),color=colour,alpha=0.8)
        plt.xlabel("Reco-Truth (cm)")
    elif Type[-1]=="X" or Type[-1]=="Y":
        MyData=MyData/36
        plt.hist(MyData[0:-1],bins=np.arange(range_start, range_end, 2),fill=False,histtype='step',range=(range_start,range_end),color=colour,alpha=0.8) 
        plt.xlabel("Reco-Truth (no. bars)")
    else:
        MyData=MyData/10
        plt.hist(MyData[0:-1],bins=60,fill=False,histtype='step',range=(range_start,range_end),color=colour,alpha=0.8)
        plt.xlabel("Reco-Truth (cm)")
    #We add parameter choice, tag, mean, median, std, and IQD to the list of handles for the legend
    patch=mpatches.Patch(color=colour,edgecolor=None,label=int_to_string(num)+" - %s"%(tag),alpha=0.8)
    MyDataOrdered=np.sort(MyData)
    mean=MyData.mean()
    plt.axvline(mean,color=colour,linestyle='--')
    patch_mean=mpatches.Patch(color=colour,fill=False,linestyle='--',label='Mean of %s: %2.2f'%(int_to_string(num),mean))
    median = MyDataOrdered[int(len(MyData)/2)]
    plt.axvline(median,color=colour,linestyle=':')
    patch_median=mpatches.Patch(color=colour,fill=False,linestyle=':',label='Median of %s: %2.2f'%(int_to_string(num),median))
    std=np.std(MyData)
    patch_std=mpatches.Patch(color='w',fill=False,label='$\sigma$ for %s: %2.2f'%(int_to_string(num),std))
    iqd=MyDataOrdered[3*int(len(MyData)/4)]-MyDataOrdered[int(len(MyData)/4)]
    patch_iqd=mpatches.Patch(color='w',fill=False,label='IQD for %s: %2.2f'%(int_to_string(num),iqd))
    handles=[]
    handles.extend([patch,patch_mean,patch_median,patch_std,patch_iqd])
    return handles

#This function plots two histograms on a single plot
#Inputs: Type - the coordinate/position we investigate e.g. UEndZ, nums - the two parameter choices to compare (expressed as 3-digit integers), range - a 2x1 array of range_start and range_end, colours - a 2x1 array of two colours
def Plotter(Type, nums, colours, range):
    #Run the DataProcessor twice to make histogram data for our two parameter choices
    DataProcessor(Type,nums[0],range[0],range[1],colours[0],"Optimum") 
    DataProcessor(Type,nums[1],range[0],range[1],colours[1],"Sub-optimum")
    #Create the joint legend for the two lines using handles, a zero line, and title using Type
    handles,labels=plt.gca().get_legend_handles_labels()
    plt.legend(handles=handles,loc='upper right',fontsize='medium',frameon=False) 
    plt.axvline(0,color='#000000',linestyle='-')
    if Type[-1]=="X" or Type[-1]=="Y":
        title = Type[:-1]+"NotZ"
    else:
        title=Type
    plt.title("%s Separation for Parameters: %s and %s"%(title,int_to_string(nums[0]),int_to_string(nums[1])))
    plt.ylabel("No. of Events")
    plt.grid(True,linestyle='--',alpha=0.4)
    #Save the plot
    plt.savefig('/uni-mainz.de/homes/bellbenj/Documents/Plots/ComparativeHistograms/%sComparativeHistogram.png'%(Type),bbox_inches='tight')
    plt.close()
    
#Info needed to construct comparative histograms
Type="End"
# indicate [optimum, suboptimum] parameter choices
nums=[417,218]
colours=['#0072B2','#D55E00']
range=[0,200]
Plotter(Type,nums,colours,range)
