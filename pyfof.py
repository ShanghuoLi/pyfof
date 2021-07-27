'''
To find velocity-coherent structures with friend-of-friend method.
'''


from __future__ import division
import os
import pickle
import numpy as np
from scipy import spatial
from astropy.io import fits
import timeit

##---------------------
def inner(x0,y0,v0,x,y,v,r,dv):
    '''
    Purpose: Judge whether the given point is in the given velocity region
    
    Inputs:
        (x0,y0,v0) -- centeral point of the area
        (x,y,v)    -- the point to be judged
        r          -- radius of the sphere
        dv         -- gradient z

    Returns:
        flag -- True or false
    '''

    flag = False
    r1 = np.sqrt((x-x0)**2+(y-y0)**2)
    dv1 = np.abs(v-v0)/r1
    if dv1<=dv:
        flag = True

    return flag


##---------------------
class fof(object):
    '''
    ---------  input ---------- \n\
    fname: fits file name.  \n\
    compt: velocity components. \n\
    I0: the minimum peak inensity of seed points. \n\
    r0 [pc]: linear separation threshold.  \n\
    dv0 [km/s/pc]: velocity gradient threshold  \n\
    advel: whether use the adaptive velocity gradient threshold, 
    if True will use the adaptive velocity gradient threshold instead of dv0. 
    You also need to provide the adaptive velocity gradient threshold in the input fits file.  \n\
    num_min: minimal number of points in a group such as number of pixels of one beam or two beam. \n\
    dist [pc]: distance \n
    --------- return ------------ \n\
    grp: the group \n\
    iso: Stage-0 isolated points \n\
    iso1: the total isolated points \n\
    '''
    

    def __init__(self, fitsname='./test.fits', compt=3, I0=0.04, dr=0.023, dv=3, num_min=150, dist=1300, advel=False, plot=True):
        self.name = fitsname
        self.compt = compt
        self.I0 = I0
        self.num_min = num_min
        self.advel = advel
        self.dist = dist

        # load data
        hdu = fits.open(self.name)
        h = hdu[0].header
        img = hdu[0].data
        s = img.shape
        th_p = h['CDELT2']

        self.dr0 = dr/self.dist*180/np.pi/th_p  # [pix]
        print('dr = {:.5f} pc = {:.5f} pix'.format(dr, self.dr0))

        if self.advel == False:
            dv0 = dv * self.dist/180*np.pi*th_p  # [km/s/pix]
            print('dv = {:.5f} km/s/pc = {:.5f} km/s/pix'.format(dv, dv0))


        #====================================================================
        #------------ First Step ----------------
        # Find the seed points friends
        #-----------------------------------------
        # change the type of img to list
        Itot = img[0:self.compt,:,:]
        vtot = img[self.compt:self.compt*2,:,:]

        if self.advel == True:
            sigtot = img[self.compt*2:self.compt*3,:,:]

            sig = []  # list of points [ra,dec,sigma]
            for i in range(s[1]):
                for j in range(s[2]):
                    for si in sigtot[:,i,j]:
                        if not np.isnan(si):
                            sig.append([si])

        tmap = []
        for i in range(s[1]):
            for j in range(s[2]):
                for v in vtot[:,i,j]:
                    if not np.isnan(v):
                        tmap.append([i,j,v])

        ## seed point candicates
        nvlsr = np.where( Itot > self.I0, vtot, np.nan)

        # change the type of img to list
        lop = []  # list of points [ra,dec,v]
        rdlop = []
        for i in range(s[1]):
            for j in range(s[2]):
                for v in nvlsr[:,i,j]:
                    if not np.isnan(v):
                        lop.append([i,j,v])
                        rdlop.append([i,j])

        
        ##------ start to run fof 
        starttime = timeit.default_timer()

        iso = []
        grp = []
        L = len(lop)
        print('Number of seed point={}'.format(L))
        start = 0
        while lop:
            nxt = []
            #if start == 0:
            #   nxt.append(lop.pop(startpoint))
            #   rdlop.pop(startpoint)
            #   start = 1
            nxt.append(lop.pop(0))
            rdlop.pop(0)

            ## run the loop to find the friends
            p = 0      
            while p < len(nxt):
                tree = spatial.KDTree(rdlop)
                ind = tree.query_ball_point([nxt[p][0], nxt[p][1]], self.dr0)

                if self.advel == True:
                    # using adaptive velocity graident
                    tind = tmap.index([ nxt[p][0], nxt[p][1], nxt[p][2] ])
                    dv0 = sig[tind][0]
                # print ('dv0', dv0)   

                if len(ind)>1:
                    tmgrp = []
                    for k in ind:
                        tmgrp.append( [ lop[k][0], lop[k][1], lop[k][2] ] )

                    ## Judging whether is a friend
                    for tlop in tmgrp:
                        if inner(nxt[p][0],nxt[p][1],nxt[p][2], tlop[0],tlop[1],tlop[2], self.dr0, dv0):
                            tempind = lop.index( tlop )
                            nxt.append(lop.pop(tempind))
                            rdlop.pop(tempind)
                            # print( 'remaining points%'%(len(lop)) )
                    print ( 'Stage-1: {:.3f}%'.format(100*(1-len(lop)/L)) )
                # print('p',p, '\n', 'len(nxt)',len(nxt))
                p += 1        
       
                ## end the inner loop if run out of the data point
                if len(lop) <= 1:
                    print('Loop ended')
                    break

            ## distinguish the group and isolate
            if len(nxt) < self.num_min:
                iso.append(nxt)
            else:
                grp.append(nxt)

            ## end the outer loop if run out of the data point
            if len(lop) <= 1:
                print('Loop ended')
                break 

        ##====================================================================
        ##------------ Second Step ----------------
        ## Find the seed point's friends
        ##-----------------------------------------

        ## change the type of img to list
        lop1 = []  # list of points [ra,dec,v]
        rdlop1 = []
        #Ilop1 = []
        for i in range(s[1]):
            for j in range(s[2]):
                for v in img[self.compt:self.compt*2,i,j]:
                    if not np.isnan(v):
                        lop1.append([i,j,v])
                        rdlop1.append([i,j])
                        #Ilop1.append([v])

        ## remove the seed data points
        for i in range(len(grp)):
            for j in range(len(grp[i])):
                tempind = lop1.index( [ grp[i][j][0], grp[i][j][1], grp[i][j][2] ] )
                lop1.pop(tempind)
                rdlop1.pop(tempind)

        L = len(lop1)
        for i in range(len(grp)):
            lop1.pop(0)
            rdlop1.pop(0)

            for j in range(len(grp[i])):
                tree = spatial.KDTree(rdlop1)
                ind = tree.query_ball_point([grp[i][j][0],grp[i][j][1]], self.dr0)

                if self.advel == True:
                    # using adaptive velocity graident
                    tind = tmap.index([ grp[i][j][0], grp[i][j][1], grp[i][j][2] ])
                    dv0 = sig[tind][0]
                # print ('dv0', dv0)   
        
                if len(ind) > 1:
                    tmgrp = []
                    for k in ind:
                        tmgrp.append( [ lop1[k][0], lop1[k][1], lop1[k][2] ] )

                    # find the friends 
                    for tlop in tmgrp:
                        if inner(grp[i][j][0],grp[i][j][1],grp[i][j][2], tlop[0],tlop[1],tlop[2], self.dr0, dv0):
                            tempind = lop1.index( tlop )
                            grp[i].append( lop1.pop(tempind) )
                            rdlop1.pop( tempind )
                            # print ('Stage-2 remaining points:: %s'%(len(lop1)) )
                    print ( 'Stage-2: {:.3f}%'.format(100*(1-len(lop1)/L)) )

                ## end the inner loop if run out of the data point
                if len(lop1) <= 1:
                    print('Loop ended')
                    break

            ## end the outer loop if run out of the data point
            if len(lop1) <= 1:
                print('Loop ended')
                break 


        ## skip this part because the remaining points unlikely have abiliity to become a group
        # #====================================================================
        # #------------ Third Step ----------------
        # #-----------------------------------------
        # # Re-make friend for isolated points 
        # #===========================================

        # iso2 = []
        # iso2 = lop1

        # rdlop2 = []
        # lop2 = []
        # for k in range(len(iso2)):
        #     lop2.append([iso2[k][0], iso2[k][1], iso2[k][2]])
        #     rdlop2.append([iso2[k][0], iso2[k][1]])


        # tmap = []
        # for i in range(s[1]):
        #     for j in range(s[2]):
        #         for v in img[4:8,i,j]:
        #             if not np.isnan(v):
        #                 tmap.append([i,j,v])

        # sig = []  # list of points [ra,dec,sigma]
        # for i in range(s[1]):
        #     for j in range(s[2]):
        #         for si in img[8:12,i,j]:
        #             if not np.isnan(si):
        #                 sig.append([si])

        # isoiso = []
        # isogrp = []
        # L = len(lop2)
        # while lop2:
        #     nxt = []
        #     nxt.append(lop2.pop(0))
        #     rdlop2.pop(0)

        #     p = 0  
        #     while p < len(nxt): 
        #         tree = spatial.KDTree(rdlop2)
        #         ind = tree.query_ball_point([nxt[p][0], nxt[p][1]], self.dr0)
        
        #         # adopt adaptive velocity
        #         tind = tmap.index([ nxt[p][0], nxt[p][1], nxt[p][2] ])
        #         dv = sig[tind][0]

        #         if len(ind)>1:
        #             tmgrp = []
        #             for k in ind:
        #                 tmgrp.append( [ lop2[k][0], lop2[k][1], lop2[k][2] ] )

        #             for tlop in tmgrp:
        #                 if inner(nxt[p][0],nxt[p][1],nxt[p][2], tlop[0],tlop[1],tlop[2], self.dr0, self.dv0):
        #                     tempind = lop2.index( tlop )
        #                     nxt.append( lop2.pop(tempind) )
        #                     rdlop2.pop( tempind )
        #                     # print ('Stage-3 remaining points:',len(lop2))
        #             print ( 'Stage-2: %.3f%%'%(100*(1-len(lop2)/L)) )

        #         p += 1

        #         ## end the inner loop if run out of the data point
        #         if len(lop) <= 1:
        #             print('Loop ended')
        #             break

        #     if len(nxt) < min_num_group:
        #         isoiso.append(nxt)
        #     else:
        #         isogrp.append(nxt)

        #     ## end the outer loop if run out of the data point
        #     if len(lop) <= 1:
        #         print('Loop ended')
        #         break 


        stoptime = timeit.default_timer()
        print('Time: {} mins = {} hrs'.format((stoptime - starttime)/60,(stoptime - starttime)/3600) )


        self.grp = grp      ## the group
        self.is0 = iso      # the isolated points in the first stage
        self.is1 = lop1     # the final isolated points


        ## output ======================================================================
        ##--------------------------
        # from datetime import datetime
        # now = datetime.now()
        # current_time = now.strftime("%Y-%m-%d_%H:%M:%S")
        
        # output = '_%s.p'%(current_time)
        # os.system('rm -rf ./variable/'+'*'+output)
        # pickle.dump(self.grp,open('./variable/grp'+output,'wb'))
        # pickle.dump(self.iso,open('./variable/iso'+output,'wb'))
        # pickle.dump(self.iso1,open('./variable/iso1'+output,'wb'))


        # ======================== analysis ==============================================
        # to screen -------------------------------------------------------------------
        l_grp = len(grp)
        n_grp = np.zeros(l_grp,dtype=np.int)  # number of points in each group
        for i in range(l_grp):
            n_grp[i] = len(grp[i])
        N_grp = sum(n_grp)  # number of points in group

        N_iso = 0
        n_iso = np.zeros(len(iso),dtype=np.int)
        for i in range(len(iso)):
            n_iso[i] = len(iso[i])
            N_iso += len(iso[i])

        loppp = []  # list of points [ra,dec,v]
        for i in range(s[1]):
            for j in range(s[2]):
                for v in img[self.compt:self.compt*2,i,j]:
                    if not np.isnan(v):
                        loppp.append([i,j,v])
        LL = len(loppp)
        # components analysis
        print('--------------------------------------------------------------')
        print('-------------- statistic data points --------------')
        print('Find {} groups.'.format(l_grp) )
        print('{} points in the groups, ratio = {}/{}={:.2f}%'.format(N_grp, N_grp, LL, N_grp/LL*100) )
        print('Stage-0 isolated points = {}, ratio = {:.2f}%'.format(N_iso, N_iso/LL*100))
        print('Stage-1 isolated points (total isolated points) = {}, ratio = {:.2f}%'.format(len(lop1),len(lop1)/LL*100))
        print('--------------------------------------------------------------')



        ###----------->>> plot 3D  <<<<---------------
        if plot==True:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            ra = []
            dec = []
            v = []
            for i in range(l_grp):
                ra.append([])
                dec.append([])
                v.append([])
                for j in range(n_grp[i]):
                    dec[i].append(grp[i][j][0])
                    ra[i].append(grp[i][j][1])
                    v[i].append(grp[i][j][2])

            color_set = ['b','k','r','cyan','fuchsia','lime','olive','darkviolet','peru', 'darkgreen', 'purple']
            fig = plt.figure()  
            ax = fig.add_subplot(111, projection='3d')

            for i in range(len(grp)):
                c = color_set[divmod(i,len(color_set))[1]]
                ax.scatter(ra[i],dec[i],v[i],marker='.',s=2, c=c)
            # plt.title('3d groups, number of groups: %d'%l_grp)
            ax.set_xlabel('Ra[pix]')
            ax.set_ylabel('Dec[pix]')
            ax.set_zlabel('v[km/s]')
            plt.show()



