#These will be useful
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#This is a function which creates a heatmap for our parameter space of extrapolation parameters
#Inputs: PerformanceData - the data to plot on the heat map (as a matrix), name - the name of the chart, Type - the coordinate/position we investigate e.g. UEndZ
def HeatmapPlotting(PerformanceData,name,Type):
    fig,ax=plt.subplots()
    #Making the plot (have to transpose just because ¯\_(ツ)_/¯ ) and set axes
    plot=ax.imshow(abs(PerformanceData.T),cmap='summer_r')
    #Set our axes - x is for ExtrapoalteDist, y for ExtrapolateLimit
    x=np.arange(1,10)
    y=np.arange(15,24)
    ax.set_xticks(range(len(x)))
    ax.set_xticklabels(x)
    ax.set_yticks(range(len(y)))
    ax.set_yticklabels(y)
    ax.set_xlabel('ExtrapolateDist')
    ax.set_ylabel('ExtrapolateLimit')
    #For labelling the minimum point - which we find below in the loop
    min_P=abs(PerformanceData[0][0])
    min_P_index="---"
    #Populate the heatmap with the information in the input matrix
    for i in range(len(x)):
        for j in range(len(y)):
            ax.text(i,j,'%.2f'%abs(PerformanceData[i][j]),ha='center',va='center',color='k',size='small')
            if abs(PerformanceData[i][j])<min_P:
                min_P=abs(PerformanceData[i][j])
                min_P_index="%d,%d"%(i+1,j+15)
    #Format legend to display min_P_index at bottom of screen
    patch=mpatches.Patch(color='w',label='Parameters minimising performance: %s'%(min_P_index))
    handles,labels=plt.gca().get_legend_handles_labels()
    handles.extend([patch])
    plt.legend(handles=handles,loc='lower center',bbox_to_anchor=(0.5,-0.3,0,0),fontsize='large',frameon=False)
    #Include a colour bar
    cbar=ax.figure.colorbar(plot,ax=ax)
    #Add label axes - either cm or bars (cm only for NotZ data)
    if Type[-1]=="X" or Type[-1]=="Y":
        unit = "no. bars"
    else:
        unit = "cm"
    cbar.ax.set_ylabel('Mean Separation (%s)'%(unit),rotation=90)
    #Add a title using the name and Type inputs
    if Type[-1]=="X" or Type[-1]=="Y":
        title = Type[:-1]+"NotZ"
    else:
        title=Type
    ax.set_title('%s at %s'%(name,title))
    #Finally, we save the plot using name and Type inputs
    fig.tight_layout()
    plt.savefig('/uni-mainz.de/homes/bellbenj/Documents/Plots/%s%s.png'%(Type,name),dpi=200,bbox_inches='tight')
    print('Created: %s%s.png'%(Type,name))
    plt.close()

#This function tells us which plots to make (histograms and heatmaps) 
#Inputs: Type - the coordinate/position we investigate e.g. UEndZ, range_start - start limit for histograms,range_end - end limit for the histograms
def Plotting(Type,range_start,range_end):
    #Create empty data sets for the heatmaps we will plot eventually (default value is an arbitratily large number)
    PerformanceData=np.ones((9,9))*9999
    PerformanceData2=np.ones((9,9))*9999
    ZThinPerformanceData=np.ones((9,9))*9999
    ZThickPerformanceData=np.ones((9,9))*9999
    ZDoubleThickPerformanceData=np.ones((9,9))*9999
    #We loop over each data set we have prepared for each choice of parameters (i - ExtrapolateDist, j - ExtrapolateLimit)
    for i in range (1,10):
        for j in range (15,24):
            #Read in the .txt files created by the Evaluation .cpp scripts. There are 81 files for each Type
            MyData=np.loadtxt("/uni-mainz.de/homes/bellbenj/Documents/Extrapolation Stuff/Histogram%sData%d%d.txt"%(Type,i,j), float)
            #Plotting the histogram - different unit for X and Y compared to Z and only have multiple lines for Z
            #For types which are EndZ we consider different stopping types (which region of steel does it stop in)
            if Type[-1]=="Z" and Type[-2]=="d":
                #Every other element of MyData is the end Z point and the rest are the reco-truth separations for that track (/10 to convert to cm)
                MyZData=MyData[1:len(MyData):2]
                MyData=MyData[0:len(MyData):2]/10
                #We split up the EndZ separations into sets for tracks which stop in: (a) thin steel region, (b) '' thick '', (c) '' double thich ''
                ZThinData=[]
                ZThickData=[]
                ZDoubleThickData=[]
                for k in range(len(MyZData)-1):
                    if MyZData[k]<14435:
                        ZThinData.append(MyData[k]) 
                    elif 14435<=MyZData[k]<17495:
                        ZThickData.append(MyData[k]) 
                    elif MyZData[k]>=17495:
                        ZDoubleThickData.append(MyData[k])
                #The heatmap data will be the means of these sets
                ZThinPerformanceData[i-1][j-15]=np.mean(ZThinData)
                ZThickPerformanceData[i-1][j-15]=np.mean(ZThickData)
                if len(ZDoubleThickData)>0:
                    ZDoubleThickPerformanceData[i-1][j-15]=np.mean(ZDoubleThickData)
                #We make a histogram with all 3 sets in it, in different colours
                counts,bins,patches=plt.hist([ZThinData,ZThickData,ZDoubleThickData],bins=60,fill=False,histtype='step',color=['#0072B2','#E69F00','#009E73'],range=(range_start,range_end))
                patch_thin=mpatches.Patch(color='#0072B2',edgecolor=None,label='Stopping in Thin Steel')
                patch_thick=mpatches.Patch(color='#E69F00',edgecolor=None,label='Stopping in Thick Steel')
                patch_doublethick=mpatches.Patch(color='#009E73',edgecolor=None,label='Stopping in Double Thick Steel')
                plt.xlabel("Reco-Truth (cm)")

                #Introduce a measure for how quickly the tail falls off
                counts_total=counts[0]+counts[1]+counts[2]
                start_mean = max(counts_total)
                for k in range(len(counts)):
                    end_mean=sum(counts_total[3*k:3*k+3])/3
                    if end_mean>=start_mean/15:
                        tail_drop=bins[3*k]
                        ratio_tail=(3*k)/len(bins)
                        plt.axvline(tail_drop,color='#D55E00',linestyle='--',label='Tail End: %d, %.2f'%(tail_drop,ratio_tail))
                        break   
            #For StartZ Types, we just plot the separation info on a histogram
            elif Type[-1]=="Z" and Type[-2]=="t":
                MyData=MyData[0:len(MyData):2]/10
                plt.hist(MyData[0:-1],bins=60,color='#0072B2',range=(range_start,range_end))
                plt.xlabel("Reco-Truth (cm)")
            #For any X or Y types, we plot the separation histograms (/36 to convert to bars)
            elif Type[-1]==("X" or "Y"):
                MyData=MyData/36
                #We have one bin per bar
                counts,bins,patches=plt.hist(MyData[0:-1],bins=range(range_start,range_end, 2),color='#0072B2',range=(range_start,range_end)) 
                plt.xlabel("Reco-Truth (no. bars)")
            #For Start or End Types, we plot histograms
            else:
                MyData=MyData/10
                counts,bins,patches=plt.hist(MyData[0:-1],bins=60,color='#0072B2',range=(range_start,range_end))
                plt.xlabel("Reco-Truth (cm)")
                #Introduce a measure for how quickly the tail falls off
                start_mean = max(counts)
                for k in range(len(bins)):
                    end_mean=sum(counts[3*k:3*k+3])/3
                    if end_mean<=start_mean/15:
                        tail_drop=bins[3*k+3]
                        ratio_tail=(3*k+3)/len(bins)
                        plt.axvline(tail_drop,color='#D55E00',linestyle='--',label='Tail End: %d, %.2f'%(tail_drop,ratio_tail))
                        break

            #Sort the data and measure length
            MyDataOrdered=np.sort(MyData)
            size = len(MyData)
            #Mean, median, quartiles, std
            mean=MyData.mean()
            median = MyDataOrdered[int(size/2)]
            first_q = MyDataOrdered[int(size/4)]
            third_q = MyDataOrdered[3*int(size/4)]
            iqd=third_q-first_q
            std=np.std(MyData)
            #Record the performance metric - try out using either median or mean as the metric
            PerformanceData[i-1][j-15]=mean 
            PerformanceData2[i-1][j-15]=median
            #Add all interesting quantities to the histogram
            plt.axvline(mean,color='k',linestyle='--',label='Mean: %2.2f'%(mean))
            plt.axvline(median,color='k',linestyle=':',label='Median: %2.2f'%(median))
            plt.axvline(0,color='#D55E00',linestyle='-',label='Zero')
            #For all but EndZ Types (for which it is too confusing) add a translucent span for IQD
            if not(Type[-1]=="Z" and Type[-2]=="d"):
                plt.axvspan(first_q,third_q,color='#CC79A7',alpha=0.5)
            #Make labels for IQD and std which can be added to the legend 
            patch_std=mpatches.Patch(color='w',fill=False,label='$\sigma$: %2.2f'%(std))
            patch_iqd=mpatches.Patch(color='#CC79A7',edgecolor=None,label='IQD: %2.2f'%(iqd),alpha=0.5)
            handles,labels=plt.gca().get_legend_handles_labels()
            if Type[-1]=="Z" and Type[-2]=="d":
                handles.extend([patch_thin,patch_thick,patch_doublethick,patch_std])
            else:
                handles.extend([patch_std,patch_iqd])
            #Label the axes and title depending on the Type and add grid and legend
            if Type[-1]=="X" or Type[-1]=="Y":
                title = Type[:-1]+"NotZ"
            else:
                title=Type
            plt.title("%s Separation for Parameters: %d,%d"%(title,i,j))
            plt.ylabel("No. of Events")
            plt.grid(True,linestyle='--',alpha=0.4)
            plt.legend(handles=handles,loc='upper right',fontsize='medium',frameon=False)
            #Save and then clear the plot before we make the next one
            plt.savefig('/uni-mainz.de/homes/bellbenj/Documents/Plots/%sHistograms/%sHistogram%d%d.png'%(Type,Type,i,j),bbox_inches='tight')
            plt.close()

    #For EndZ only, we make a heatmap for stopping in: thin, thick, and double-thick layers
    if Type[-1]=="Z" and Type[-2]=="d":
        HeatmapPlotting(ZThinPerformanceData,"ZThinPerformance",Type)
        HeatmapPlotting(ZThickPerformanceData,"ZThickPerformance",Type)
        HeatmapPlotting(ZDoubleThickPerformanceData,"ZDoubleThickPerformance",Type)
    #For all types, we make heatmaps of mean and median
    HeatmapPlotting(PerformanceData,"Performance",Type)
    HeatmapPlotting(PerformanceData2,"Performance2",Type)

#We make plots for various different Types (uncomment as needed)
    #These are the Types from after the 3d track matching:
Plotting("Start",0,100)
Plotting("End",0,300)
    #These are the Types from the LineCandidate info:
#Plotting("UStartX",-20,20)
#Plotting("UStartZ",-50,50)
Plotting("UEndX",-50,50)
Plotting("UEndZ",-200,250)
#Plotting("VStartX",-10,10)
#Plotting("VStartZ",-50,50)
#Plotting("VEndX",-50,50)
Plotting("VEndZ",-200,250)
#Plotting("XStartY",-10,10)
#Plotting("XStartZ",-50,50)
#Plotting("XEndY",-50,50)
Plotting("XEndZ",-200,250)
