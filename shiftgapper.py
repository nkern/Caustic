import numpy as np
def shiftgapper(data,sort=True):
    npbin = 25
    gap_prev = 2000.0
    nbins = np.int(np.ceil(data[:,0].size/(npbin*1.0)))
    origsize = data[:,0].shape[0]
    data = data[np.argsort(data[:,0])] #sort by r to ready for binning
    print 'NBINS FOR GAPPER = ', nbins
    for i in range(nbins):
        print 'BEGINNING BIN:',str(i)
        databin = data[npbin*i:npbin*(i+1)]
        print 'R BETWEEN', str(databin[:,0][0]),'and',str(databin[:,0][-1])
        print 'DATA SIZE IN',databin[:,0].size
        datanew = None
        nsize = databin[:,0].size
        datasize = nsize-1
        if nsize > 5:
            while nsize - datasize > 0 and datasize >= 5:
                print '    ITERATING'
                nsize = databin[:,0].size
                databinsort = databin[np.argsort(databin[:,1])] #sort by v
                f = (databinsort[:,1])[databinsort[:,1].size-np.int(np.ceil(databinsort[:,1].size/4.0))]-(databinsort[:,1])[np.int(np.ceil(databinsort[:,1].size/4.0))]
                gap = f/(1.349)
                print '    GAP SIZE', str(gap)
                if gap < 500.0: break
                #if gap >= 2.0*gap_prev:
                #    gap = gap_prev
                #    print '   Altered gap = %.2f'%(gap)
                databelow = databinsort[databinsort[:,1]<=0]
                gapbelow =databelow[:,1][1:]-databelow[:,1][:-1]
                dataabove = databinsort[databinsort[:,1]>0]
                gapabove = dataabove[:,1][1:]-dataabove[:,1][:-1]
                try:
                    if np.max(gapbelow) >= gap: vgapbelow = np.where(gapbelow >= gap)[0][-1]
                    else: vgapbelow = -1
                    print 'MAX BELOW GAP',np.max(gapbelow)
                    try: 
                        datanew = np.append(datanew,databelow[vgapbelow+1:],axis=0)
                    except:
                        datanew = databelow[vgapbelow+1:]
                except ValueError:
                    pass
                try:
                    if np.max(gapabove) >= gap: vgapabove = np.where(gapabove >= gap)[0][0]
                    else: vgapabove = 99999999
                    print 'MAX ABOVE GAP',np.max(gapabove)
                    try: 
                        datanew = np.append(datanew,dataabove[:vgapabove+1],axis=0)
                    except:
                        datanew = dataabove[:vgapabove+1]
                except ValueError:
                    pass
                databin = datanew
                datasize = datanew[:,0].size
                datanew = None
            print 'DATA SIZE OUT', databin[:,0].size
            if gap >= 500.0:
                gap_prev = gap
            else:
                gap_prev = 500.0
        try:
            datafinal = np.append(datafinal,databin,axis=0)
        except:
            datafinal = databin
    print 'GALAXIES CUT =',str(origsize-datafinal[:,0].size)
    if sort == True:
        datafinal = datafinal[np.argsort(datafinal[:,2])]
    return datafinal
