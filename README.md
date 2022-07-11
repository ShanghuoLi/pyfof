# pyfof
python-based friend-of-friend algorithm to identify the filament in the position-position-velocity space. 

You can simply run the pyfof.py in the ipython

> from pyfof.py import fof 
> 
> res = fof(fitsname='./test_for_pyfof.fits', compt=3, I0=0.04, dr=0.023, dv=3, num_min=150, dist=1300, advel=False, plot=True) 
> 
Then, you can manually save the results to the pick files.
> 
> output = 'pyfof_results.p' 
>
> pickle.dump(res.grp,open('grp_'+output,'wb')) 
> 
> pickle.dump(res.iso,open('iso_'+output,'wb')) 
> 
> pickle.dump(res.iso1,open('iso1_'+output,'wb')) 


fname: fits file name. 

compt: velocity components. 

I0: the minimum peak inensity of seed points. 

r0 [pc]: linear separation threshold.  

dv0 [km/s/pc]: velocity gradient threshold.

advel: whether use the adaptive velocity gradient threshold, if True will use the adaptive velocity gradient threshold instead of dv0. 
You also need to provide the adaptive velocity gradient instead of sigma in the input fits file.  

num_min: minimal number of points in a group such as number of pixels of one beam or two beam. 

dist [pc]: distance.

plot: plot the 3d scatter. 


