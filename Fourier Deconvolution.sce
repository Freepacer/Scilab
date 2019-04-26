clear;
clf(0); clf(1); clf(2);clf(3);clf(4);
tic();

function y=fgls(x,a)
    y = a(1)*(exp( -4*log(2)*(x-a(3))^2/(a(2))^2 ));
endfunction

function e=myfun(a,xm,ym,wm)
    e = wm.*(fgls(xm,a) - ym)
endfunction

    x1=-9.736;
    dx=x1-0.9096;
    
    filename='scan1 -9.736-510.555.txt'; // enter file name here
    data=fscanfMat(filename);
    [n col]=size(data);
    disp('filename: '+filename)

    xm=(data(:,2)+dx)*8.06573; ym=data(:,3);

    [ymax,iymax]=max(ym);
    
    //use 60 points for fitting; 40 starting from peak and 20 last points
    xmf=zeros(60,1);
    xmf(1:40)=xm(iymax:iymax+39);
    xmf(41:60)=xm(n-19:n);
    ymf(1:40)=ym(iymax:iymax+39);
    ymf(41:60)=ym(n-19:n);
    wmf=ymf./ymf;

    // separate fit to get peak center
    xmfp=zeros(23,1);
    xmfp=xm(iymax-8:iymax+14);
    ymfp=ym(iymax-8:iymax+14);
    wmfp=ymfp./ymfp;

    a0= [ymax; 50; 0]    // initial guess
    [fpopt,apopt,gpopt]=leastsq(list(myfun,xmfp,ymfp,wmfp),a0);
    [fopt,aopt,gopt]=leastsq(list(myfun,xmf,ymf,wmf),a0);
    
    yfit=fgls(xm,aopt);
    yfitp=fgls(xm,apopt);
    
    disp(aopt)

    disp('peak center:')
    disp(apopt(3))

    scf(0)
    plot2d(xm,ym,-1);
    plot2d(xm,yfit,2);
    plot2d(xm,yfitp,5);
    
    xm=xm-apopt(3);

    //scf(1)
    //plot2d(xm,ym,-1);
    //plot2d(xmf-apopt(4),ymf,-2); // show points used in fit
    //plot2d(xm,yfit,2); // show fit function
    //plot2d(xm,yfitp-aopt(5),5); // show fit function for peak center
    
    //scf(2)
    
    //plot2d(xm,ym-yfit);

    ftym=fft(ym);
    ftyfit=fft(yfitp);
    
    
    decon=ifft(log(ftym./ftyfit).*ftyfit);
    
    output=[ym decon];
    
    scf(3);
    plot2d(xm,ym,2);
    plot2d(xm,decon,5);
    //disp('Lorentzian Percentage:');
    //LP=apopt(4)*100;
    //disp(LP);
    disp(apopt(1));
    disp(apopt(2));
    disp(apopt(3));
    
    scf(4);
    //plot2d(xm,ftyfit,2);
    //plot2d(xm,ftym,5);
    plot2d(xm,ftym./ftyfit,2);
    
    outfilename=filename+'_out.txt';
    ofd=mopen(outfilename,'w');
    mfprintf(ofd,'%f\t%f\n',output);
    mclose(ofd);
toc();
