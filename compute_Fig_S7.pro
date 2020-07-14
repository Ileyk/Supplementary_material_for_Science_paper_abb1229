pro compute_Fig_S7


  a_au = 15.
  M1 = 1.02
  M2 = 0.51
  vwind_terminal = 15.
  
  
  !p.thick=3
  !p.charthick=3

  ;set_plot,'ps'
  ;device,/color,file='Fig_S7.eps',/encapsulated,yoff=0,/portrait, xoff=0,xsize=15,ysize=10, set_font='Courier',/tt_font

  grav = 6.67259d-8 ; cm3 g-1 s-2
  M_SUN    =1.989d+33 ;gram
  AU       =1.496d13 ; cm
  R_sun    =6.96d10 ; cm
  ysec     =31557600d0 ; sec

  a = a_au * au ; binary separation in cm

  ;;;;;; PLOT ;;;;
  gridding=100.
  n_elem = floor(2*a/au)*gridding
  x=(findgen(n_elem)+1)/gridding ; in au
  xh=where(x ge 1.) ; restrict above 1au
  x=x(xh); in au; restrict above 1au
  ; compute escape velocity from M1 for plotting
  vesc = sqrt(2*grav*M1*M_sun/(x*au)) / 1d5 ; escape velocity in km/s
  max_vesc = max(vesc)

  ymin = 0.
  ymax=max_vesc
  plot, x/a_au, vesc, xtitle='distance / a', ytitle='velocity [km/s]',xra=[0.,2.],/xs,yra=[ymin, ymax],/ys
  
  ;;;;;;;

  ; compute orbital period
  orbital_period = sqrt(4 * !Pi^2 * a^3./ (grav * (M1*M_sun+M2*M_sun))) ; in sec
  print, 'orbital period is ', orbital_period/ysec,' year'

  
  ; compute displacement of star1 w.r.t. center of mass
  r1_G = (M2/(M1+M2))*a ; in cm
  ; compute displacement of star2 w.r.t. center of mass
  r2_G = (M1/(M1+M2))*a ; in cm

  print, 'distance of star1 to COM is ', r1_G/au, ' au'
  print, 'distance of star2 to COM is ', r2_G/au, ' au'
   ver, r1_G/a, line=1,color=fsc_color('deep pink')
  xyouts, r1_G/a,max_vesc/1.05,'CoM',alignment=-0.1,color=fsc_color('deep pink'),charsize=0.7

  ; compute place of L1 w.r.t. star2 (less massive one); if M2<<M1
  L1 = (M2/3.*M1)^(1./3.)*a ; in cm
  print, 'distance of L1 w.r.t. star2 is ', L1/au,' au if M2/M1<<1'
  



  ; calculate the radius of the Roche Lobe of each star 
  ; in units of the orbital separation
  q = M2/M1
  r_RL2 = 0.49 * q^(2./3.)/(0.6*q^(2./3.)+alog(1.+q^(1./3.))) ; in units of a
  q = M1/M2
  r_RL1 = 0.49 * q^(2./3.)/(0.6*q^(2./3.)+alog(1.+q^(1./3.))) ; in units of a
  
  R_RL2_cm = R_RL2*a
  R_RL1_cm = R_RL1*a

  print, 'radius of Roche Lobe star1 is ', r_RL1, ' a'
  print, 'radius of Roche Lobe star2 is ', r_RL2, ' a'
  arrow, 0,39,+R_RL1,39., thick=3, /data,/solid,hthick=1
  xyouts,0.1, 40, 'R!DR,1!N',charsize=0.7,alignment=-0.5
  arrow, 1.-R_RL2,39,1.+R_RL2,39., thick=3, color=fsc_color('gray'), /data,/solid,hthick=1
  arrow, 1.+R_RL2,39,1.-R_RL2,39., thick=3, color=fsc_color('gray'), /data,/solid,hthick=1
  xyouts,1., 40, 'R!DR,2!N',charsize=0.7,color=fsc_color('gray'),alignment=-0.5
  ver,1.


  ;;;;;;;;;;;;;
  ; SPEEDS    ;
  ;;;;;;;;;;;;;
  ; for circular orbit: compute orbital speed of star2 around star1 (hence at binary separation)
  orbital_speed = 2 * !pi*a/orbital_period ; in cm/s
  print, 'orbital speed at distance of binary separation is ', orbital_speed/1d5,' km/s'
  plotsym,3,1,/fill ; filled star,  default size
  oplot,a_au[0:0]/a_au,orbital_speed[0:0]/1d5,psym=8
  xyouts,a_au/a_au, orbital_speed/1d5, 'v!Dorb!N(m!Dcomp!N w.r.t. M!D*!N)',alignment=-0.1,charsize=0.7

  ; for circular orbit: compute orbital speed of star1 circling around center-of-mass
  orbital_speed = 2 * !pi*r1_G/orbital_period ; in cm/s
  print, 'orbital speed of star1 w.r.t. COM is ', orbital_speed/1d5,' km/s'
  oplot,x[0:0]/a_au,orbital_speed[0:0]/1d5,psym=8
  xyouts,x[0]/a_au, orbital_speed/1d5, 'v!Dorb!N(M!D*!N w.r.t. CoM)',alignment=-0.1,charsize=0.7



  ; for circular orbit: compute orbital speed of star2 circling around center-of-mass
  orbital_speed = 2 * !pi*r2_G/orbital_period ; in cm/s
  print, 'orbital speed of star2 w.r.t. COM is ', orbital_speed/1d5,' km/s'
  oplot,a_au[0:0]/a_au,orbital_speed[0:0]/1d5,psym=8
  xyouts,a_au/a_au, orbital_speed/1d5, 'v!Dorb!N(m!Dcomp!N w.r.t. CoM)',alignment=-0.1,charsize=0.7

  ; for circular orbit: compute orbital speed at roche lobe of star1 ; assuming same orbital period
  orbital_speed = 2 * !pi*r_RL1*a/orbital_period ; in cm/s
  print, 'orbital speed at roche lobe of star1 is ', orbital_speed/1d5,' km/s'

  ; for circular orbit: compute orbital speed at roche lobe of star2 ; assuming same orbital period
  orbital_speed = 2 * !pi*r_RL2*a/orbital_period ; in cm/s
  print, 'orbital speed at roche lobe of star2 is ', orbital_speed/1d5,' km/s'

  ; escape velocity of star1 at roche lobe radius
  escape_velocity = sqrt(2*grav*M1*M_sun/(R_RL1*a))
  print, 'escape velocity of star1 at roche lobe radius is ', escape_velocity/1d5, ' km/s'

  
  beta_c=0.5 ; for C-rich AGB wind
  beta_o=5. ; for O-rich AGB wind
  
  v_infty = vwind_terminal ; in km/s
  v_0 = 2. ; in km/s
  R_star = 267.*R_sun ; in cm from Chen et al. 2019
  R_dust = 2.75 *R_star; in cm
  vwind_o=fltarr(n_elements(x))
  vwind_c=fltarr(n_elements(x))
  xt=where(x*au ge R_dust)
  vwind_o(xt) = v_0 + (v_infty - v_0)*(1.- R_dust/(x(xt)*au))^beta_o ; in cm/s
  vwind_c(xt) = v_0 + (v_infty - v_0)*(1.- R_dust/(x(xt)*au))^beta_c ; in cm/s
  xt=where(x*au lt R_dust)
  vwind_o(xt) = v_0 
  vwind_c(xt) = v_0 
  ver, R_dust/a,line=3
  polyfill, [R_star/a, R_star/a,0,0],[ymin, ymax, ymax, ymin], /data, orientation = 45
  ver, R_star/a


  oplot,x/a_au,vwind_o,color=fsc_color('red')
  oplot,x/a_au,vwind_c,color=fsc_color('blue')
  xyouts,1.6, max(vwind_o), 'v!Dwind!N(O-rich)',alignment=-0.1,charsize=0.8,color=fsc_color('red')
  xyouts,1.6, max(vwind_c)*1.03, 'v!Dwind!N(C-rich)',alignment=-0.1,charsize=0.8,color=fsc_color('blue')


  ;device,/close
  ;set_plot,'x'


end