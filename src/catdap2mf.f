cc      program catdap
c
c      categorical data analysis program package 02  (catdap 02 )
c
c                                                     july  1, 1982
c                                              revised oct. 1, 1988
c                                              revised jan. 14, 1993
c                                              revised Jun. 24, 1999
c                                              revised by makio ishiguro  aug. 6, 2015
c
c
c           designed by        yosiyuki sakamoto
c           coded by   kouichi katsura and yosiyuki sakamoto
c
c                   the institute of statistical mathematics
c                   4-6-7 minami-azabu, minato-ku, tokyo 106
c
c
c          this program simultaneously searches for the best subset and
c       the best categorization of explanatory variables which have the
c       most effective information on a specific categorical variable.
c       this program is applicable to any case provided that the
c       response variable is categorical.
c          for the datails of "how to operate catdap-02", see appendix
c       of "categorical data analysis by aic" (sakamoto (1991)).
c          note that this program requires a work file with device
c       reference number '1'.
c
c     ----------------------------------------------------------------
c
c     < how to write a data descriptive file >
c
c     examples of the data descriptive file are given in sample input
c     (a) to (c) in appendix of sakamoto(1991).
c
c     row 1:(in cols 1-72)
c       title of the data descriptive file.
c
c     row 2:(in cols 1-72)
c       title of main data set.
c
c     row 3:(in cols 1-72)
c       blank should be inserted between entries.
c       <entry 1>: sample size
c       <entry 2>: number of variables
c       <entry 3>: number of variables having code values to be
c                  converted
c       <entry 4>: number of variables used for the specification
c                  for neglect of a part of the main data, if any
c       <entry 5>: type of the main data
c                  =0: integer
c                  =1: otherwise
c       <entry 6>: =1: to perform only the analysis of two-way tables
c                  =0: otherwise
c       <entry 7>: number of variables to be retained for the analysis
c                  of multidimensional tables
c       <entry 8>: number of subsets of variables if additional
c                  analyses for any subsets are to be made
c       <entry 9>: number of two-way tables in the output <3>
c                  ('0' is regarded as 10. )
c      <entry 10>: =1: to obtain the output <4>
c                  =0: otherwise
c
c     row 4-5:(in cols 1-72)
c       read format of the main data.  always use 2 rows.
c
c     row 6 and beyond:
c       the following input <entry 1-9> is required every variable.
c       <entry 1>:  (in cols 1-3)
c                variable number
c       <entry 2>:  (in cols 6-25)
c                title of the variable
c       <entry 3>:  (in cols 26-28)
c                minimum of the code values
c       <entry 4>:  (in cols 29-31)
c                maximum of the code values
c       <entry 5>:  (in cols 32-34)
c                the way of pooling of categories of the variable;
c            =0: equally spaced pooling
c            =1: unequally spaced pooling
c            =2: no pooling
c       <entry 6>:  (in cols 35-37)
c            =1: if the variable is considered as the response variable
c            =0: otherwise
c       <entry 7>:  (in cols 38-40)
c            =1: if any code values of the variable should be
c                converted.
c               (also, specify new code values corresponding to the
c                original code values 1, 2, ... , 19 in the following
c                row.   for example, see lines 8 and 10 in sample
c                input(a) of appendix of sakamoto(1991).)
c            =0: otherwise
c       <entry 8>:  (in cols 41-43)
c            =1: if the variable pertains to neglect of a part of the
c                main data.  for example, see line 10 in sample
c                input(c) of appendix of sakamoto(1991).
c
c                this program discards a sample in which the variable
c                takes a value within one of several intervals (or takes
c                a code value) to be specified by (a) (or by (b)) below.
c
c               (a) if and only if this <entry 8> is 1 and <entry 9> is
c                not 0.0, specify the number of intervals to be ignored
c                in the next row, and the left and right end values in
c                the row after next.  for example, see lines 11 and 12
c                in sample input (c) of appendix of sakamoto(1991).
c               (b) if this <entry 8> is 1, and if <entry 9> is 0.0 or i
c                <entry 5> of row 3 is o, specify the number of code
c                values and code values to be ignored in the next row.
c
c            =0: otherwise
c
c       <entry 9>:(in cols 44-46)
c            the accuracy of measurement
c            (set 0.0 if the variable is an integer.)
c
c     the last row(s):(in cols 1-72)
c        the number of variables and those variable numbers for each
c        of which the additional analysis is to be made, if <entry 8>
c        of row 3 is not 0
c
c     ----------------------------------------------------------------
c
c     < possible output >
c         <1> data descriptive file
c         <2> list of single explanatory variables (arranged in
c             ascending order of aic)
c         <3> the corresponding two-way tables (each having an optimal
c             categorization)
c         <4> graphical representation of the two-way tables
c         <5> aic's of the models with k explanatory variables
c             (k=1,2,...)
c                     and
c             summary of subsets of explanatory variables
c         <6> contingency table constructed by the best subset of
c             explanatory variables
c         <7> the output of the additional analysis
c
c     ----------------------------------------------------------------
c
cc      implicit real*8(a-h,o-z)
cc      integer *2 ia
cc      integer alim,recode
cc      character* 4 fmt
cc      dimension ia(900000),a(55000),fmt(100)
ccc     dimension idt(8)
ccc     dimension datet(20),dataf1(20),dataf2(20)
ccc     character*24 datet,dataf1,dataf2,dataf3,dataf4
cc      character*24 dataf1,dataf2,dataf3,dataf4
cc      character*10 datet(3)
cc      common ialim,imax,alim,jmax,imaxx,jmaxx
ccc     data datet/20*'    '/
ccc
ccc     data input
ccc
cc      ialim=900000
cc      alim=55000
cc      open (4,file='catdap02.out')
ccc     open (2,file='om94.dsc')
cc      open (1,file='catdap.temp',form='unformatted')
cc      write(4,2007)
cc      write(6,*) 'input dsc filename'
cc      read(5,*) dataf3
cc      write(6,*) 'input main data filename'
cc      read(5,*) dataf4
cc      open (2,file=dataf3)
ccc     read(5,1004,end=200) (dataf1(i),i=1,10)
ccc     read(5,1004,end=200) (dataf2(i),i=1,10)
cc      read(2,1005,end=200) dataf1
cc      read(2,1005,end=200) dataf2
ccc     write(4,2008) (dataf1(i),i=1,10),(dataf2(i),i=1,10)
cc      write(4,2008) dataf1,dataf2
ccc     call fdate(datet)
ccc     call date(datet)
ccc     call date_and_time(datet(1),datet(2),datet(3),idt)
ccc     write(4,2009) (datet(i),i=1,10)
cc      write(4,2009) (datet(i),i=1,2)
cc      write(4,2010)
ccc     write(4,2014) (dataf1(i),i=1,10)
ccc     write(4,2014) (dataf2(i),i=1,10)
cc      write(4,1005) dataf1
cc      write(4,1005) dataf2
cc      in=10
ccc     open(in,file='om94r1')
cc      open(in,file=dataf4)
cc      read(2,*,end=200) nsamp,n,recode,iskip1,it,ipart,nov,icl,ikp,izu
cc      write(4,2011) nsamp,n,recode,iskip1,it,ipart,nov,icl,ikp,izu
cc      nfm=2
cc      nfmm=nfm*18
cc      read(2,1004,end=200) (fmt(i),i=1,nfmm)
cc      mmm=nfm
cc      do 151 j=1,mmm
cc      m1=1+18*(j-1)
cc      m2=18*j
cc      write(4,2014) (fmt(i),i=m1,m2)
cc  151 continue
cc      if(in.eq.0) in=5
cc      if(nov.eq.0) nov=n
cc      l=1
cc      samp=nsamp
cc      n22=2*n+1
cc      n24=4*n
cc      mmm=(2*n-1)/20+1
cc      do 110 j=1,mmm
cc      m1=n22+20*(j-1)
cc      m2=n22+20*j-1
cc      if(m2.gt.n24) m2=n24
cc  110 continue
cc      n2=0
cc      n2a=n22
cc      n22=0
cc      do 10 i=1,n
cc      n23=2*i+n2a-1
cc      inn=i+n
cc      ia(i)=ia(n23-1)
cc      ia(inn)=ia(n23)
cc      n1=ia(inn)-ia(i)+1
cc      n2=max(n2,n1)
cc      n22=n22+n1
cc   10 continue
cc      if(it.ne.1) go to 120
cc      mmm=(n-1)/20+1
cc      do 130 j=1,mmm
cc      m1=1+20*(j-1)
cc      m2=20*j
cc      if(m2.gt.n) m2=n
cc  130 continue
cc  120 continue
cc      i1=n+1
cc      i2=i1+n
cc      i0=2*n
cc      i0n=i0+n
cc      i01=i0+1
cc      mmm=(n-1)/20+1
cc      do 140 j=1,mmm
cc      m1=i01+20*(j-1)
cc      m2=i01+20*j-1
cc      if(m2.gt.i0n) m2=i0n
cc  140 continue
cc      i3=i0+n+1
cc      i0=i3-1
cc      if(recode.eq.0) go to 15
cc      do 20 k=1,recode
cc      i01=i0+1
cc      i02=i0+20
cc      do 11 i=i01,i02
cc      if(ia(i).eq.0) go to 12
cc   11 continue
cc      i=i02+1
cc   12 continue
cc      i0=i0+20
cc   20 continue
cc   15 continue
cc      i4=i3+20*recode
cc      i5=i4+l
cc      nfmm=nfm*20
cc      mmm=nfm
cc      do 150 j=1,mmm
cc      m1=1+20*(j-1)
cc      m2=20*j
cc  150 continue
cc      i0=i5-1
cc      do 40 k=1,n
cc      i01=i0+1
cc      i0=i0+20
cc   40 continue
cc      i6=i5+20*n
cc      i0=i6-1
cc      if(iskip1.eq.0) go to 50
cc      do 45 k=1,iskip1
cc      i01=i0+1
cc      i02=i0+20
cc      do 41 i=i01,i02
cc      if(ia(i).eq.0) go to 42
cc   41 continue
cc      i=i02+1
cc   42 continue
cc      i0=i0+20
cc   45 continue
cc      if(it.eq.1) go to 50
cc      i0=i7-1
cc      j0=n
cc      do 56 k=1,iskip1
cc      j00=j0
cc      i01=i0+1
cc      i02=i0+2
cc      i0=i02
cc      j0=j00+20
cc   56 continue
cc      j1=n+1
cc   50 continue
cc      i7=i6+20*iskip1
cc      i8=i7
cc      if(it.eq.1) i8=i7+2*iskip1
cc      if(icl.eq.0) go to 90
cc      i01=i8
cc      i02=i01+9
cc      do 100 k=1,icl
cc      do 91 i=i01,i02
cc      if(ia(i).eq.0) go to 92
cc   91 continue
cc      i=i02+1
cc   92 continue
cc      i01=i02+1
cc      i02=i01+9
cc  100 continue
cc   90 continue
cc      i9=i8+icl*10
ccc
ccc     storage allocation
ccc
cc      j1=1
cc      if(it.eq.1) j1=1+n
cc      call input(n,recode,iskip1,it,icl,nov,
cc     1           ia(1),ia(i1),ia(i2),ia(i3),ia(i4),ia(i5),ia(i6),ia(i7),
cc     2           ia(i8),a(1),a(j1),in)
cc
cc      write(4,2017) nsamp
cc      write(4,2018) n
cc      write(4,2019) recode
cc      write(4,2022) iskip1
cc      write(4,2023) it
cc      write(4,2024) ipart
cc      write(4,2025) nov
cc      write(4,2029) icl
cc      write(4,2033) ikp
cc      write(4,2034) izu
cc      write(4,2016)
cc      write(4,2026)
cc      write(4,2016)
cc      write(4,2027)
cc      write(4,2028)
cc      write(4,2030)
cc      write(4,2031)
cc      write(4,2032)
cc      iaf=ia(i4)
cc      iaff=ia(i4)+n
cc      n11=ia(iaff)-ia(iaf)+1
cc      if(ia(iaff).eq.0.and.ia(iaf).eq.0) then
cc      write(4,2036)
cc      stop 10
cc 2036 format(' catdap can be applied to data sets in which response' 
cc     &       ' variables are categorical.')
cc      end if
cc      n22=0
cc      do 35 i=1,n
cc      inn=i+n
cc      n1=ia(inn)-ia(i)+1
cc      n2=max(n2,n1)
cc      n22=n22+n1
cc   35 continue
cc      n33=100
cc      if(nsamp/2+1.lt.n33) n33=nsamp/2+2
cc      n33=max(n33,n2)
cc      if(it.eq.1) n22=n33*n
cc      if(it.eq.1) go to 30
cc      n33=n2
cc   30 continue
cc      i10=i9+n
cc      i11=i10+n*n33
cc      i12=i11+n
cc      i13=i12+n
cc      i14=i13+n33
cc      i15=i14+n11*n22
cc      i16=i15+n11*n22
cc      i17=i16+10*n33
cc      i18=i17+n11
cc      i19=i18+n33
cc      i20=i19+n
cc      i21=i20+n33
cc      i22=i21+n33
cc      i23=i22+n33
cc      i24=i23+10
cc      i25=i24+10
cc      i26=i25+n
cc      i27=i26+n
cc      imaxx=i27
cc      if(it.eq.0) go to 70
cc      j2=j1+20*iskip1
cc      j3=j2+n
cc      j4=j3+n*n33
cc      j5=j4+n33
cc      j6=j5+n11*n33
cc      j7=j6+n11*n33
cc      j8=j7+n11*n33
cc      j9=j8+n33
cc      j10=j9+n11
cc      j11=j10+n
cc      j12=j11+n
cc      j13=j12+n
cc      j14=j13+n
cc      j15=j14+n
cc      j16=j15+n33
cc      j17=j16+n33
cc      j18=j17+n33
cc      go to 80
cc   70 j1=1
cc      j2=j1+20*iskip1
cc      j3=j2+n
cc      j4=j3
cc      j5=j4+n33
cc      j6=j5+n11*n33
cc      j7=j6+n11*n33
cc      j8=j7+n11*n33
cc      j9=j8+n33
cc      j10=j9+n11
cc      j11=j10
cc      j12=j11
cc      j13=j12
cc      j14=j13
cc      j15=j14
cc      j16=j15+n33
cc      j17=j16+n33
cc      j18=j17+n33
cc   80 continue
cc      jmaxx=j18
cc      iskip1=iskip1+1
cc      recode=recode+1
cc      icl=icl+1
cc      nnn=ialim-i16
cc      jnn=alim-j4
cc      inn=ialim-i27
cc      jjn=alim-j18
cc      if(imaxx.gt.ialim.or.jmaxx.gt.alim) go to 60
cc      imax=i14
cc      jmax=j4
cc      call cat2(nsamp,n,l,recode,in,iskip1,it,fmt,ipart,nov,izu,icl,
cc     1          ikp,ia(1),ia(i1),ia(i2),ia(i3),ia(i4),ia(i5),
cc     2          ia(i6),ia(i7),ia(i8),ia(i9),ia(i10),ia(i11),ia(i12),
cc     3          ia(i13),ia(i14),ia(i15),ia(i16),ia(i17),ia(i18),ia(i19),
cc     4          ia(i20),ia(i21),ia(i22),ia(i23),ia(i24),ia(i25),ia(i26),
cc     5          a(1),a(j1),a(j2),a(j3),a(j4),
cc     6          a(j5),a(j6),a(j7),a(j8),a(j9),a(j10),a(j11),a(j12),
cc     7          a(j13),a(j14),a(j15),a(j16),a(j17),samp,n11,n22,
cc     8          n33,ia(i27),a(j18),nnn,jnn,inn,jjn)
cc      write(4,2006) imaxx,jmaxx
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop
cc   60 write(4,2005) ialim,imaxx,alim,jmaxx
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
cc  200 write(4,2035)
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
cc 1004 format(20a4)
cc 1005 format(a24)
cc 2001 format(' ')
cc 2005 format(' ia or a dimension over ',4i10)
cc 2006 format(2x,'used memory'//
cc     1       i10,2x,'words  ( for ''ia'' )'/i10,2x,
cc     2       'words  ( for ''a'' )')
cc 2007 format(72('*')/'*',70x,'*'/
cc     1          '*  ***      *     *****   ****      *     ',
cc     2                  '****          ***     ***    *'/
cc     3          '* *   *    * *      *      *  *    * *    ',
cc     4                  '*   *        *   *   *   *   *'/
cc     5          '* *       *   *     *      *  *   *   *   ',
cc     6                  '*   *        *   *       *   *'/
cc     7          '* *       *   *     *      *  *   *   *   ',
cc     8                  '****         *   *    ***    *'/
cc     9          '* *       *****     *      *  *   *****   ',
cc     a                  '*     *****  *   *   *       *'/
cc     b          '* *   *   *   *     *      *  *   *   *   ',
cc     c                  '*            *   *   *       *'/
cc     d          '*  ***    *   *     *     ****    *   *   ',
cc     e                  '*             ***    *****   *'/
cc     f          '*',70x,'*',/'*',70x,'*'/
cc     g          '*    categorical   data   analysis   progra',
cc     h                  'm   package   02            *'/
cc     i          '*',70x,'*'/72('*'))
ccc2008 format(/44x,'revised version oct 1, 1988'/
ccc    1         22x,'designed by y.sakamoto'/
ccc    2         13x,'programmed by k.katsura and y.sakamoto'/
ccc    3         12x,'the institute of statistical mathematics'/
ccc    4         12x,'4-6-7 minami-azabu, minatoku, tokyo 106 '///
ccc    5         'title of data descriptive file : ',10a4/
ccc    6         'title of data set              : ',10a4)
cc 2008 format(/44x,'revised version oct 1, 1988'/
cc     1         22x,'designed by y.sakamoto'/
cc     2         13x,'programmed by k.katsura and y.sakamoto'/
cc     3         12x,'the institute of statistical mathematics'/
cc     4         12x,'4-6-7 minami-azabu, minatoku, tokyo 106 '///
cc     5         'title of data descriptive file : ',a24/
cc     6         'title of data set              : ',a24)
ccc2009 format(  'date                           : ',10a4)
cc 2009 format(  'date                           : ',a8,a10)
cc 2010 format(/'<1>'/'+',23('-'),'+'/'i data descriptive file i'/
cc     &       '+',23('-'),'+'/72('-'))
cc 2011 format(14i5)
cc 2012 format(7f10.4)
cc 2013 format(10a1)
cc 2014 format(18a4)
cc 2015 format(1h+,101x,'!')
cc 2016 format(72('-'))
cc 2017 format('<note>'/'row  1;(in cols 1-72)'/
cc     1                '  title of the data descriptive file.'/
cc     2                'row  2;(in cols 1-72)'/
cc     3                '  title of the main data set.'/
cc     4        72('-')/'row  3;(in cols 1-72)'/
cc     5                '  blank should be inserted between entries.'/
cc     6                '  <entry 1>: sample size =',i5)
cc 2018 format('  <entry 2>: number of variables =',i3)
cc 2019 format('  <entry 3>: number of variables having code values to be'
cc     &      /'             converted =',i3)
cc 2020 format('  number of data to be skipped to required data =',i3)
cc 2021 format('  input data device number =',i3)
cc 2022 format('  <entry 4>: number of variables used for the specificatio
cc     &n'/    '             for neglect of a part of the main data, ifan
cc     &y'/    '             =',i3)
cc 2023 format('  <entry 5>: type of main data( integer=0, otherwise=1 )
cc     &=',i3)
cc 2024 format('  <entry 6>: to perform only the analysis of two-way table
cc     &s(=1)'/'             or otherwise(=0) =',i3)
cc 2025 format('  <entry 7>: number of variables to be retained for the an
cc     &alysis'/'             of multidimensional tables =',i3)
cc 2029 format('  <entry 8>: number of subsets of variables if additional
cc     &analyses'/'             for any subsets are to be made =',i3)
cc 2033 format('  <entry 9>: number of two-way tables in the output <3> '
cc     1      /'             (''0'' is regarded as 10. ) =',i3)
cc 2034 format('  <entry 10>: =1: to obtain the output <4>'/
cc     1       '              =0: otherwise     =',i3)
cc 2026 format('row 4-5;(in cols 1-72)'/
cc     &       '  read format of the main data. always use 2 rows.')
cc 2027 format('row 6 and beyond:'/
cc     1       '  the fllowing input <entry 1-9> is required every variabl
cc     2e.'/
cc     3       '  <enrry 1>:  (in cols 1-3)'/
cc     4       '           variable number'/
cc     5       '  <entry 2>:   (in cols 6-25)'/
cc     6       '           title of the variable '/
cc     7       '  <entry 3>:   (in cols 26-28)'/
cc     8       '           minimum of the code values'/
cc     9       '  <entry 4>:   (in cols 29-31)'/
cc     a       '           maximum of the code values'/
cc     b       '  <entry 5>:   (in cols 32-34)'/
cc     c       '           the way of pooling categories of the variable;'
cc     d      /'       =0: equally spaced pooling'/
cc     e       '       =1: unequally spaced pooling'/
cc     f       '       =2: no pooling')
cc 2028 format('  <entry 6>:   (in cols 35-37)'/
cc     1       '       =1: if the variable is considered as the response v
cc     2ariable'/
cc     3       '       =0: otherwise'/
cc     4       '  <entry 7>:   (in cols 38-40)'/
cc     5       '       =1: if any code values of the variable should be'/
cc     6       '           converted.'/
cc     7       '          (also, specify new code values corresponding to
cc     8the'/
cc     9       '           original code values 1,2,...,19 in the followin
cc     ag'/
cc     b       '           row.  for example, see lines 8 and 10 in sample
cc     c'     /'           inpit(a) of appendix of sakamoto(1991).)'/
cc     d       '       =0: otherwise')
cc 2030 format('  <entry 8>:   (in cols 41-43)'/
cc     1       '       =1: if the variable pertains to the neglect of apa
cc     2rt'/
cc     3       '           of the main data.  for example, see line 10 in
cc     4sample'/
cc     5       '           input(c) of appendix of sakamoto(1991).)'//
cc     6       '           this program discards a sample in which the var
cc     7iable'/
cc     8       '           takes a value within one of several intervals (
cc     9or takes'/
cc     a       '           a code value) to be specified by (a) (or by (b)
cc     b) below.'//
cc     c       '          (a) if and only if this <entry 8> is 1 and <entr
cc     dy 9> is')
cc 2031 format('           not 0.0, specify the number of intervals to be
cc     1ignored'/
cc     2       '           in next row, and the left and right end values
cc     3in'/
cc     4       '           the row after next.  for example, see lines 11
cc     5and 12'/
cc     6       '           in sample input(c) of appendix of sakamoto(1991
cc     7).'/
cc     8       '          (b) if this <entry 8> is 1, and if <entry 9> is
cc     90.0 or if'/
cc     a       '           <entry 5> of row 3 is 0, specify the number of
cc     bcode'/
cc     c       '           values and coce values to be ignored in the nex
cc     dt row.'/
cc     e       '       =0: otherwise')
cc 2032 format('  <entry 9>:   (in cols 44-46)'/
cc     1       '       the accuracy of measurement'/
cc     2       '       (set 0.0 if the variable is an intrger.)'/
cc     3       72('-')/
cc     4       'the last row(s):   (in cols 1-72)'/
cc     5       '   the number of variables and those variable numbers of e
cc     6ach'/
cc     7       '   of which the additional analysis is to be made, if <ent
cc     8ry 8>'/
cc     9       '   of row 3 is not 0')
cc 2035 format(' end of file error')
cc      end
cc      subroutine cat2(nsamp,n,l,recode,in,iskip1,it,fmt,ipart,nov,izu,
cc     1                icl,ikp,item1,item2,ity,iconv,face,
cc     2                title,iskip,isk,icls,ite,iaa,item,idata,
cc     3                totalr,ia,iab,lcy,totalc,ttrr,kp,nc,ncc,knn,lj,lt,
cc     4                idat,noo,xx,sk,dx,ab,tttr,a1,a,pa,ptr,ptc,
cc     5                da,bmin,bmax,am1,ai,aa,ax,totr,samp,n11,n22,
cc     6                n33,iw,w,nw,jnn,inn,jjn)
       subroutine catdap2mf(nsamp,n,l,recode,iskip1,it,ipart,nov,icl,
     1              item1,item2,ity,iconv,face,iskip,isk,sk,icls,xx,
     2              ida,da,iab,totalc,ttrr,ab,iaa,pa,idata,ite,dx,
     3              aaam,caa,icaa,nnia,lk77,morder,iby,ibc,pbc,
     4              aic1,iabse,baseaic,n11,n33,ikr,jkr,ikkk,eps01,ier)
c
c     this subroutine is substantially the main program of catdap-02.
c     this searches for maice within the possible ways of categorization
c     of a single explanatory variable in case 'ity(k)=0'.
c
      INCLUDE 'catdap_f.h'
c
      implicit real*8(a-h,o-z)
cc      integer*2 item,ia,idata,totalr,face,ite,title,iab,iaa,lcy,iw,isk,
cc     1          totalc,ttrr,nc,ncc,ity,kp,knn,lj,lt,idat,item1,item2,
cc     2          iconv,iskip,noo,icls,bl,mx,st
      integer totalr,face,totalc,ttrr,caa
      integer recode
cc      character*4 fmt
c
cc      dimension item1(n),item2(n),item(n),face(l),title(20,n),fmt(100),
cc     1          idata(nw),tttr(jnn),ia(n11,n22),iskip(20,iskip1),noo(n),
cc     2          totalr(n33),iconv(20,recode),totr(n33),ttrr(n33),dx(n),
cc     3          a(n11,n33),iab(n11,n22),iaa(n,n33),lcy(10,n33),knn(n33),
cc     4          totalc(n11),a1(n11,n33),pa(n11,n33),ptr(n33),ptc(n11),
cc     5          da(n),bmin(n),bmax(n),xx(n),am1(n),ai(n),kp(n),aa(n33),
cc     6          ab(n,n33),nc(n33),ncc(n33),isk(2,iskip1),sk(20,iskip1),
cc     7          lj(10),ax(n33),lt(10),ite(n),ity(n),idat(n),iw(inn),
cc     8          w(jjn),icls(10,icl)
      dimension item1(n),item2(n),ity(n),iconv(20,recode),face(l),
     1           iskip(20,iskip1),isk(2,iskip1),sk(20,iskip1),
     2           icls(max(10,nov+1),icl),xx(n),ida(n,nsamp),da(n,nsamp),
     3           iab(n,n11,n33),totalc(n11),ttrr(n,n33),ab(n,n33),
     4           iaa(n,n33),pa(n,n11,n33),ptr(n,n33),ptc(n11),idata(n),
     5           ite(n),dx(n),aaam(ikr),caa(ikr,jkr),
     6           icaa(ikr),nnia(ikr),morder(ikr),ier(3),
     7           iby(nov,ikkk,icl+1),ibc(n11,ikkk,icl+1),
     8           pbc(n11,ikkk,icl+1),aic1(icl+1),
     9           iabse(nov,n33,icl+1)
      dimension item(n),totalr(n33),ia(n,n11,n33),tttr(n33),a1(n11,n33),
     1           lcy(10,n33),kp(n),nc(n33),ncc(n33),knn(n33),lj(10),
     2           lt(10),idat(n),a(n11,n33),bmin(n),bmax(n),am1(n),ai(n),
     3           aa(n33),ax(n33),totr(n33)
c
cc      character*4    fm1(20),fm2(20),fm3(20)
cc      character*1    bl,mx,st
cc      common ialim,imax,jalim,jmax,imaxx,jmaxx
cc      data bl,mx,st/' ','-','*'/
cc      data (fm1(i),i=1,20)
cc     &        /'    ','    ',',2x,','    ','f6.1',',1x,','1h(,','i4, ',
cc     &         '1h))',11*'    '/
cc      data (fm2(i),i=1,20)
cc     &        /'   2','   3','   4','   5','   6','   7','   8','   9',
cc     &         '  10','  11',10*'    '/
cc      data (fm3(i),i=1,20)
cc     &        /'(i5 ','    ','(5ht','otal',16*'    '/
c
c
c      continuous response variable
c
      magic = 0
      ires = face(1)
      if (xx(ires) .ne. 0) then
         magic = ires
         imagic = item2(ires)
      end if
c
c     initialization
c
c <<<
      ikp=n
      ig=face(1)
      if(ig.eq.0) return

      iab = 0
      totalc = 0
      ttrr = 0.0d0
      ab = 0.0d0
      iaa = 0
      pa = 0.0d0
      idata = 0
      ite = 0
      dx = 0.0d0
      aaam = 0.0d0
      caa = 0.0d0
      icaa = 0
      nnia = 0
      morder = 0
      iby = 0
      ibc = 0
      pbc = 0
      aic1 = 0
      iabse = 0
      ier = 0
c >>>
      nx1=n
cc      iskip1=iskip1-1
cc      recode=recode-1
      do 40 i=1,n11
cc      do 40 j=1,n22
cc      ia(i,j)=0
cc   40 iab(i,j)=0
      do 39 j=1,n33
      do 38 k=1,n
      ia(k,i,j)=0      
   38 iab(k,i,j)=0
   39 continue
   40 continue
      if(it.eq.0) go to 80
c
c     choice of the initial class interval for continuous data in
c     the case of 'it=1'
c
cc      ig=face(1)
      nj=nsamp/2
      nns=0
      do 90 j=1,nsamp
cc      read(in,fmt) (da(i),i=1,n)
      if(iskip1.eq.0) go to 85
      do 86 i=1,iskip1
      iskip2=isk(1,i)
      if(xx(iskip2).eq.0.0) go to 88
      iskip3=isk(2,i)
      do 87 ii=1,iskip3
cc      if(da(iskip2).gt.sk(2*ii-1,i).and.da(iskip2).lt.sk(2*ii,i))
      if(da(iskip2,j).gt.sk(2*ii-1,i).and.da(iskip2,j).lt.sk(2*ii,i))
     1 go to 90
   87 continue
      go to 86
   88 continue
      iskip3=iskip(2,i)
cc      idata1=da(iskip2)+0.1
      idata1=da(iskip2,j)
      do 89 ii=1,iskip3
      if(idata1.eq.iskip(ii+2,i)) go to 90
   89 continue
   86 continue
   85 continue
      nns=nns+1
cc      write(1) (da(i),i=1,n)
c <<<
      if(j.gt.nns) then
         do 84 i=1,n
            da(i,nns)=da(i,j)
   84    continue
      end if
c >>>
      ii=0
      do 100 i=1,n
      if(xx(i).eq.0) go to 100
      if(nns.ne.1) go to 110
cc      bmin(i)=da(i)
cc      bmax(i)=da(i)
      bmin(i)=da(i,j)
      bmax(i)=da(i,j)
      go to 100
  110 continue
cc      if(bmin(i).ge.da(i)) bmin(i)=da(i)
cc      if(bmax(i).le.da(i)) bmax(i)=da(i)
      if(bmin(i).ge.da(i,j)) bmin(i)=da(i,j)
      if(bmax(i).le.da(i,j)) bmax(i)=da(i,j)
  100 continue
   90 continue
      nsamp=nns
      samp=nns
c <<<
      nj=nsamp/2
c >>>
      do 120 i=1,n
c-----    modified by M.I.
      if(i .eq. magic) then
           ity (i) = 2 
           bmagic = bmin(magic)
           xmagic = (bmax(magic)-bmin(magic))/item2(magic)
           do 115 kk=0,item2(magic)
             ab(magic, kk+1) = bmin(magic)+xmagic*kk
  115      continue
           xx(i) = 0.d0
           goto 120
      endif
c-----
      if(xx(i).eq.0) go to 120
      mm=(bmax(i)-bmin(i))/xx(i)+1
      mm=min(mm,50,nj)
      m=((bmax(i)-bmin(i))/xx(i)  )/mm+1
      ai(i)=xx(i)*m
      i1=bmin(i)/ai(i)
cc------------------------------ cat1
cc      am1(i)=(i1+0.5)*ai(i)
cc------------------------------
      am1(i)=i1*ai(i)+0.5*xx(i)
      if(am1(i).ge.bmin(i)) am1(i)=am1(i)-ai(i)
      if(abs(am1(i)-bmin(i)).lt.xx(i)/1.d5) am1(i)=am1(i)-ai(i)
      am2=am1(i)
      ab(i,1)=am2
      do 130 ii=1,100
      am2=am2+ai(i)
      ab(i,ii+1)=am2
      if(am2.ge.bmax(i)) go to 140
  130 continue
  140 kp(i)=ii
  120 continue
cc      rewind 1
   80 continue
c
c     main data input  (in the case of 'it=0')
c
      nns=0
c
c     construction of two-way tables
c
      do 150 k=1,nsamp
      if(it.eq.0) go to 160
c
c     discretizing continuous data into categorical data in the case
c     of 'it=1'
c
cc      read(1) (da(j),j=1,n)
      do 170 j=1,n
c-----    modified by M.I.
         if(j .eq. magic) then
           iv = (da(j,k)-bmagic)/xmagic + 1
           if(iv .gt. imagic) iv = imagic
           da(j,k) = iv
         endif
c-----
c     if(xx(j).eq.0) go to 180
      if(item1(j).ne.0.or.item2(j).ne.0) go to 180
      kpp=kp(j)
      item(j)=kp(j)
      ite(j)=kp(j)
      am11=am1(j)
      aii=ai(j)
      do 190 ii1=1,kpp
      am11=am11+aii
      if(ig.eq.j) am11=ab(ig,ii1+1)
cc      if(da(j).gt.am11) go to 190
cc      if(abs(da(j)-am11).lt.xx(j)/1.d5) go to 190
cc      idata(j)=ii1
      if(da(j,k).gt.am11) go to 190
      if(abs(da(j,k)-am11).lt.xx(j)/1.d5) go to 190
      ida(j,k)=ii1
      go to 200
  190 continue
  200 continue
      go to 170
  180 continue
cc      idata(j)=da(j)+0.01
      ida(j,k)=da(j,k)+0.01
      if(k.ne.1) go to 170
      item(j)=item2(j)-item1(j)+1
      iii=1
      if(item(j).le.50) go to 206
      iii=(item(j)-1)/50
      item(j)=(item(j)-1)/iii+1
  206 continue
      ite(j)=item(j)
      item3=item(j)
      do 205 i=1,item3
      iaa(j,i)=i+item1(j)-1
  205 continue
  170 continue
c
c     code value conversion
c
      if(recode.eq.0) go to 245
      do 255 i=1,n
      do 265 ik=1,recode
      if(iconv(1,ik).ne.i) go to 265
cc      ik1=idata(i)+1
cc      if(idata(i).eq.0) ik1=11
cc      idata(i)=iconv(ik1,ik)
      ik1=ida(i,k)+1
      if(ida(i,k).eq.0) ik1=11
      ida(i,k)=iconv(ik1,ik)
      go to 255
  265 continue
  255 continue
  245 continue
c
c     off-code check:  (off-code values are converted into the max of
c     code value  'item2(i)' )
c
      do 275 i=1,n
c     if(xx(i).ne.0.) go to 275
      if(item1(i).eq.0.and.item2(i).eq.0) go to 275
c     if(idata(i).gt.item2(i).or.idata(i).lt.item1(i)) idata(i)=item2(i)
cc      if(idata(i).gt.item2(i).or.idata(i).lt.item1(i)) go to 940
      if(ida(i,k).gt.item2(i).or.ida(i,k).lt.item1(i)) go to 940
c     if(idata(i).gt.item2(i).or.idata(i).lt.item1(i)) then
c     write(4,2035) i,idata(i)
c     idata(i)=item2(i)
c     end if
cc      idata(i)=idata(i)-item1(i)+1
      ida(i,k)=ida(i,k)-item1(i)+1
  275 continue
      do 286 i=1,n
      if(item2(i)-item1(i)+1.le.50) go to 286
      iii=(item2(i)-item1(i))/50+1
cc      idata(i)=(idata(i)-1)/iii+1
      ida(i,k)=(ida(i,k)-1)/iii+1
  286 continue
      nns=nns+1
      go to 210
  160 continue
cc      read(in,fmt) (idata(i),i=1,n)
      if(k.ne.1) go to 225
      do 220 j=1,n
      item(j)=item2(j)-item1(j)+1
      iii=1
      if(item(j).le.50) go to 226
      iii=(item(j)-1)/50
      item(j)=(item(j)-1)/iii+1
  226 continue
      ite(j)=item(j)
      item3=item(j)
      do 230 i=1,item3
      iaa(j,i)=i*iii+item1(j)-1
  230 continue
  220 continue
  225 continue
c
c     deletion of a part of the main data in the case of 'it=0'
c
      if(iskip1.eq.0) go to 212
      do 215 is1=1,iskip1
      iskip2=iskip(1,is1)
      iskip3=iskip(2,is1)
      do 217 is2=1,iskip3
cc      if(idata(iskip2).eq.iskip(is2+2,is1)) go to 150
      if(ida(iskip2,k).eq.iskip(is2+2,is1)) go to 150
  217 continue
  215 continue
  212 continue
      nns=nns+1
c
c     code value conversion
c
      if(recode.eq.0) go to 240
      do 250 i=1,n
      do 260 ik=1,recode
      if(iconv(1,ik).ne.i) go to 260
cc      ik1=idata(i)+1
cc      if(idata(i).eq.0) ik1=11
cc      idata(i)=iconv(ik1,ik)
      ik1=ida(i,k)+1
      if(ida(i,k).eq.0) ik1=11
      ida(i,k)=iconv(ik1,ik)
      go to 250
  260 continue
  250 continue
  240 continue
c
c     off-code check:  (off-code values are converted into the max of
c     code value  'item2(i)' )
c
      do 270 i=1,n
cc      if(idata(i).gt.item2(i).or.idata(i).lt.item1(i)) go to 940
      if(ida(i,k).gt.item2(i).or.ida(i,k).lt.item1(i)) go to 940
c     if(idata(i).gt.item2(i).or.idata(i).lt.item1(i)) idata(i)=item2(i)
cc      idata(i)=idata(i)-item1(i)+1
      ida(i,k)=ida(i,k)-item1(i)+1
  270 continue
cc      write(1) (idata(i),i=1,n)
      do 285 i=1,n
      if(item2(i)-item1(i)+1.le.50) go to 285
      iii=(item2(i)-item1(i))/50+1
cc      idata(i)=(idata(i)-1)/iii+1
      ida(i,k)=(ida(i,k)-1)/iii+1
  285 continue
  210 continue
c
c     construction of two-way tables
c
cc      ii=0
      do 290 i2=1,l
      i3=face(i2)
      if(item1(i3).eq.0.and.item2(i3).eq.0) go to 930
cc      i=idata(i3)+ii
cc      i=ida(i3,k)+ii
cc      jj=0
      i=ida(i3,k)
      do 300 j3=1,n
cc      j2=j3
cc      j=idata(j2)+jj
cc       ia(i,j)=ia(i,j)+1
cc  300 jj=jj+item(j2)
cc  290 ii=ii+item(i3)
      j=ida(j3,k)
      ia(j3,i,j)=ia(j3,i,j)+1
  300 continue
  290 continue
  150 continue
c
c     computation of the marginal frequencies of the two-way
c     tables, and that of the corresponding percentages
c
      nsamp=nns
      samp=nns
      i1=1
cc      i3=face(1)
cc      i2=item(i3)
      do 310 ik=1,l
c <<<
      i3=face(ik)
      i2=item(i3)
c >>>
cc      j1=1
cc      j2=item(1)
cc      jjj=0
c-----    modified by M.I.
      ibaseflag = 0
c-----
      do 320 k21=1,n
      k=k21
cc <<<
      j1=1
      j2=item(k)
c >>>
      if(k.eq.i3) go to 330
      itype=ity(k)
cc      jj=0
      do 340 j=j1,j2
cc      jj=jj+1
cc  340 totalr(jj)=0
  340 totalr(j)=0
cc      ii=0
      do 350 i=i1,i2
cc      jj=0
cc      ii=ii+1
cc      totalc(ii)=0
      totalc(i)=0
      do 360 j=j1,j2
cc      jj=jj+1
cc      a(ii,jj)=ia(i,j)
cc      a1(ii,jj)=ia(i,j)
cc      totalr(jj)=totalr(jj)+ia(i,j)
cc      totalc(ii)=totalc(ii)+ia(i,j)
      a(i,j)=ia(k,i,j)
      a1(i,j)=ia(k,i,j)
      totalr(j)=totalr(j)+ia(k,i,j)
      totalc(i)=totalc(i)+ia(k,i,j)
  360 continue
  350 continue
c-----   modified by M.I.
      if( ibaseflag .ne. 1) then
         ibaseflag = 1
         baseaic = 0.d0
         sum = 0.d0
         do 365 i=i1,i2
            if( totalc(i) .gt. 0.d0) then
             tem = totalc(i) 
             baseaic = baseaic + tem * log (tem) - 1.d0
             sum = sum + totalc( i ) 
            endif
  365    continue
         baseaic =  baseaic - sum * log(sum) + 1.d0
         baseaic = -2.d0 * baseaic
         if( magic .ne.0) then
            shift = 2.d0 * sum * LOG( xmagic )
            baseaic = baseaic + shift
         endif
      endif
c-----
cc      ii=0
      itot=0
      do 370 i=i1,i2
cc      ii=ii+1
cc      jj=0
      do 380 j=j1,j2
cc      jj=jj+1
cc      if(totalr(jj).eq.0 ) pa(ii,jj)=0.
cc      if(totalr(jj).eq.0 ) go to 380
cc      pa(ii,jj)=(ia(i,j)*100.)/totalr(jj)
      if(totalr(j).eq.0 ) pa(k,i,j)=0.
      if(totalr(j).eq.0 ) go to 380
      pa(k,i,j)=(ia(k,i,j)*100.)/totalr(j)
  380 continue
cc      ptc(ii)=(totalc(ii)*100.)/samp
cc      if(totalc(ii).ne.0) itot=itot+1
      ptc(i)=(totalc(i)*100.)/samp
      if(totalc(i).ne.0) itot=itot+1
  370 continue
      if(itot.lt.2) then
cc      write(4,2037)
      ier(1)=2037
cc      stop 10
      return
 2037 format('caution:'/
     &       '    the program catdap cannot deal with data sets where',
     &       ' the number of'/
     &       'non-zero frequency categories of the response variables',
     &       ' is less than 2.')
      end if
cc      jj=0
      do 390 j=j1,j2
cc      jj=jj+1
cc      if(totalr(jj).eq.0 ) ptr(jj)=0.
cc      if(totalr(jj).eq.0 ) go to 390
cc      ptr(jj)=(totalr(jj)*100.)/totalr(jj)
      if(totalr(j).eq.0 ) ptr(k,j)=0.
      if(totalr(j).eq.0 ) go to 390
      ptr(k,j)=(totalr(j)*100.)/totalr(j)
  390 continue
      pt=(nsamp/samp)*100.
      expo=exp(-1.)
      lq=1
      if(j2-j1.eq.0) go to 330
      if(itype.lt.0) dx(k21)=1.d10
      if(itype.lt.0) go to 330
      if(itype.eq.1.and.j2-j1.gt.1) go to 430
c     if(itype.ne.2) go to 435
      if(itype.eq.0) go to 435
c
c     searching for maice within possible ways of categorization of
c     a single explanatory variable; (in the case of 'ity(k)=0')
c
      kk1=1
      kk5=j2-j1+1
      do 445 j=1,kk5
      lcy(kk1,j)=1
  445 continue
cc      jj=0
      tsmp=0.0
      do 405 j=j1,j2
cc      jj=jj+1
cc      if(totalr(jj).eq.0) go to 405
cc      tttr(jj)=0.
      if(totalr(j).eq.0) go to 405
      tttr(j)=0.
cc      ii=0
      do 400 i=i1,i2
cc      ii=ii+1
cc      if(totalc(ii).eq.0) go to 400
cc      if(a1(ii,jj).le.0.) tsmp=tsmp+expo
cc      if(a1(ii,jj).le.0.) tttr(jj)=tttr(jj)+expo
      if(totalc(i).eq.0) go to 400
      if(a1(i,j).le.0.) tsmp=tsmp+expo
      if(a1(i,j).le.0.) tttr(j)=tttr(j)+expo
  400 continue
  405 continue
      al=0.
cc      do 411 i=1,ii
      do 411 i=i1,i2
      tt=0.0
      if(totalc(i).eq.0) go to 411
cc      do 415 j=1,jj
cc      if(totalr(j).eq.0) go to 415
      do 415 j=j1,j2
      if(totalr(j).eq.0) go to 415
      if(a1(i,j).le.0.0) tt=tt+expo
  415 continue
      tt=tt+totalc(i)
cc      do 410 j=1,jj
cc      if(totalr(j).eq.0) go to 410
      do 410 j=1,j2
      if(totalr(j).eq.0) go to 410
      aaa=a1(i,j)
      if(aaa.le.0.0) aaa=expo
         al0=al
      al=al+aaa*log(aaa/(tt*(totalr(j)+tttr(j)))*(samp+tsmp))
  410 continue
  411 continue
cc      ii0=ii
cc      do 412 i=1,ii
      ii0=i2
      do 412 i=1,i2
  412 if(totalc(i).eq.0) ii0=ii0-1
cc      jj0=jj
cc      do 413 j=1,jj
      jj0=j2
      do 413 j=1,j2
  413 if(totalr(j).eq.0) jj0=jj0-1
      inn=ii0*jj0-jj0-ii0+1
      aic=-2.*(al-inn)
      dx(k21)=aic
      go to 790
  435 kkm=k
      ng=0
      kg=0
      n1=0
cc      nj=jj-1
cc      ii0=ii
cc      do 414 i=1,ii
      nj=j2-1
      ii0=i2
      do 414 i=i1,i2
  414 if(totalc(i).eq.0) ii0=ii0-1
  440 continue
      n1=n1+1
cc      do 450 jj3=1,jj
      do 450 jj3=1,j2
      totr(jj3)=0.
cc      do 450 ii1=1,ii
      do 450 ii1=1,i2
  450 a1(ii1,jj3)=0.
      jj1=0
c <<<
      jj=1
      kkj=1
c >>>
cc      do 460 j=1,jj
      do 460 j=1,j2
c <<<
      jj=j
c >>>
      nc(j)=n1
      do 460 kj=1,n1
c <<<
      kkj=kj
c >>>
      jj1=jj1+1
      totr(j)=totalr(jj1)+totr(j)
cc      do 470 i=1,ii
      do 470 i=1,i2
      a1(i,j)=a1(i,j)+a(i,jj1)
  470 continue
cc      if(jj1.eq.jj) go to 480
      if(jj1.eq.j2) go to 480
  460 continue
  480 continue
cc      nc(j)=kj
      nc(jj)=kkj
      kk=j
      al=0.
      tsmp=0.0
      do 490 j=1,kk
      if(totr(j).eq.0.) go to 490
      tttr(j)=0.
cc      do 500 i=1,ii
      do 500 i=1,i2
      if(totalc(i).eq.0) go to 500
      if(a1(i,j).le.0.0) tsmp=tsmp+expo
      if(a1(i,j).le.0.) tttr(j)=tttr(j)+expo
  500 continue
  490 continue
cc      do 510 i=1,ii
      do 510 i=1,i2
      tt=0.0
      if(totalc(i).eq.0) go to 510
      do 515 j=1,kk
      if(totr(j).eq.0.) go to 515
      if(a1(i,j).le.0.0) tt=tt+expo
  515 continue
      tt=tt+totalc(i)
      do 520 j=1,kk
      if(totr(j).eq.0.) go to 520
      aaa=a1(i,j)
      if(aaa.le.0.0) aaa=expo
      al=al+aaa*log(aaa/(tt*(totr(j)+tttr(j)))*(samp+tsmp))
  520 continue
  510 continue
      kk0=kk
      do 511 j=1,kk
  511 if(totr(j).eq.0.) kk0=kk0-1
      inn=ii0*kk0-kk0-ii0+1
      aic=-2.*(al-inn)
      aicmm=aic
      kk5=kk
      do 530 i=1,kk
      ncc(i)=nc(i)
  530 continue
      if(n1.eq.1) go to 540
      n10=n1-1
  550 continue
cc      do 560 jj3=1,jj
      do 560 jj3=1,j2
      totr(jj3)=0.
cc       do 560 ii1=1,ii
      do 560 ii1=1,i2
  560 a1(ii1,jj3)=0.
      jj1=0
      if(n10.eq.0) go to 540
      nc(1)=n10
      do 570 kj=1,n10
      totr(1)=totr(1)+totalr(kj)
cc      do 570 i=1,ii
      do 570 i=1,i2
      a1(i,1)=a1(i,1)+a(i,kj)
  570 continue
      jj1=n10
cc      do 580 j=2,jj
      do 580 j=2,j2
      nc(j)=n1
      do 580 kj=1,n1
      jj1=jj1+1
      totr(j)=totr(j)+totalr(jj1)
cc      do 590 i=1,ii
      do 590 i=1,i2
      a1(i,j)=a1(i,j)+a(i,jj1)
  590 continue
      if(jj1.eq.j2) go to 600
  580 continue
  600 continue
      kk=j
      nc(j)=kj
      al=0.
      tsmp=0.0
      do 611 j=1,kk
      if(totr(j).eq.0.) go to 611
      tttr(j)=0.
cc      do 610 i=1,ii
      do 610 i=1,i2
      if(totalc(i).eq.0) go to 610
      if(a1(i,j).le.0.0) tsmp=tsmp+expo
      if(a1(i,j).le.0.) tttr(j)=tttr(j)+expo
  610 continue
  611 continue
cc      do 620 i=1,ii
      do 620 i=1,i2
      if(totalc(i).eq.0) go to 620
      tt=0.0
      do 625 j=1,kk
      if(totr(j).eq.0.) go to 625
      if(a1(i,j).le.0.0) tt=tt+expo
  625 continue
      tt=tt+totalc(i)
      do 630 j=1,kk
      if(totr(j).eq.0.) go to 630
      aaa=a1(i,j)
      if(aaa.le.0.0) aaa=expo
      al=al+aaa*log(aaa/(tt*(totr(j)+tttr(j)))*(samp+tsmp))
  630 continue
  620 continue
      kk0=kk
      do 621 j=1,kk
  621 if(totr(j).eq.0.) kk0=kk0-1
      inn=ii0*kk0-kk0-ii0+1
      aic=-2.*(al-inn)
ccxx      if(aicmm.lt.aic) go to 640
         daic=dabs(aicmm-aic)
         if(daic.ne.0) daic=daic/max(dabs(aicmm), dabs(aic))
         if(aicmm.lt.aic .and. daic.gt.eps01) go to 640

      aicmm=aic
      kk5=kk
      do 650 i=1,kk
      ncc(i)=nc(i)
  650 continue
  640 continue
      n10=n10-1
      go to 550
  540 continue
      aa(n1)=aicmm
      knn(n1)=kk5
      if(kg.ne.kk5) go to 660
      if(ax(ng).ge.aicmm) ax(ng)=aicmm
      if(ax(ng).ge.aicmm) go to 670
      go to 680
  660 ng=ng+1
      ax(ng)=aicmm
      kg=kk5
      if(ng.le.10) go to 690
      do 700 i=2,10
      i7=i-1
      ljj=lj(i)
      k1=knn(ljj)
      do 710 j=1,k1
      lcy(i7,j)=lcy(i,j)
  710 continue
      lj(i7)=lj(i)
      lt(i7)=lt(i)
  700 continue
  670 continue
      if(ng.le.10) go to 690
      k1=knn(n1)
      do 720 j=1,k1
      lcy(10,j)=ncc(j)
  720 continue
      lt(10)=ng
      lj(10)=n1
      go to 680
  690 continue
      do 730 i=1,kk5
      lcy(ng,i)=ncc(i)
  730 continue
      lj(ng)=n1
      lt(ng)=ng
  680 continue
      if(ng.le.5) go to 740
      if(ax(ng).gt.ax(ng-5).and.ax(ng-1).gt.ax(ng-5).and.ax(ng-2).
     1     gt.ax(ng-5).and.ax(ng-3).gt.ax(ng-5).and.ax(ng-4).gt.
     2    ax(ng-5)) go to 750
  740 continue
      if(nj.eq.n1) go to 750
      go to 440
  750 continue
      am=10.**10
c <<<
      kk1=1
c >>>
      do 760 i=1,ng
ccxx      if(am.lt.ax(i)) go to 760
         damx=dabs(am-ax(i))
         if(damx.ne.0) damx=damx/max(dabs(am), dabs(ax(i)))
      if(am.lt.ax(i) .and. damx.gt.eps01) go to 760

      am=ax(i)
      kk1=i
  760 continue
      dx(k21)=am
      do 770 j=1,10
      kk2=lt(j)
      if(kk2.eq.kk1) kk3=lj(j)
      if(kk2.eq.kk1) kk5=knn(kk3)
      if(kk2.eq.kk1) go to 780
  770 continue
  780 kk1=j
      k=kkm
      go to 790
c
c     searching for maice within possible ways of categorization of
c     a single explanatory variable; (in the case of 'ity(k)=1')
c
  430 continue
cc      i11=1+n33
cc      i12=i11+n33
cc      imax1=imaxx+i12
cc      jmax1=jmaxx
cc      if(imax1.gt.ialim) go to 10
cc      imax0=imaxx
cc      jmax0=jmaxx
cc      imax=imax1
cc      if(imax1.gt.imaxx) imaxx=imax1
cc      call ac1p(a,ii,totalr,jj,iw(1),nc,ncc,lcy,lj,knn,iw(i11),a1,
cc     1          totr,tttr,aa,am,kk5,kk2,iw(i12),w(1),n11,n33,totalc)
      call ac1p(a,i2,totalr,j2,nc,ncc,lcy,lj,knn,a1,
     1          totr,tttr,aa,am,kk5,kk2,n11,n33,totalc,eps01)
cc      imaxx=imax0
cc      jmaxx=jmax0
      kk1=kk2
      dx(k21)=am
  790 continue
c
c     printing out the two-way table with maice
c
      do 801 i=i1,i2
cc      jj=j1-1
cc      j11=jjj
      jj=0
      do 801 j=1,kk5
      kj5=lcy(kk1,j)
cc      j11=j11+1
      do 801 ij1=1,kj5
      jj=jj+1
cc  801 iab(i,j11)=iab(i,j11)+ia(i,jj)
  801 iab(k21,i,j)=iab(k21,i,j)+ia(k21,i,jj)
      jity=0
cc      j11=jjj
      do 803 j=1,kk5
      iitt=0
      do 802 i=i1,i2
cc  802 iitt=iitt+iab(i,j11+j)
  802 iitt=iitt+iab(k21,i,j)
      if(iitt.eq.0) go to 803
      jity=jity+1
  803 continue
      if(jity.le.1) ity(k21)=-1
      if(jity.le.1) dx(k21)=1.d9
      if(jity.le.1) go to 870
cc      noo(k21)=jjj
      if(it.eq.0.or.xx(k).eq.0.) go to 850
      am2=am1(k)
      do 860 j=1,kk5
      kj5=lcy(kk1,j)
      am2=am2+ai(k)*kj5
      ab(k,j+1)=am2
  860 continue
      ab(k,1)=bmin(k)
      ab(k,kk5+1)=bmax(k)
      go to 870
  850 continue
      iii=(item2(k)-item1(k))/50+1
      kj5=0
      do 880 j=1,kk5
      kj5=lcy(kk1,j)+kj5
      iaa(k,j )=kj5*iii+item1(k)-1
      if(j.eq.kk5) iaa(k,j)=item2(k)
  880 continue
  870 continue
ccc      jjj=jjj+kk5
      ite(k21)=kk5
  330 continue
cc      j1=j2+1
cc      if(k21.eq.n) go to 320
cc      j2=j2+item(k21+1)
  320 continue
      do 331 i=1,n
      idat(i)=0
  331 continue
      idata(1)=i3
      idat(i3)=1
      ikn=item(i3)
      nm1=n-1
c <<<
      nw1=0
c >>>
      do 332 i=1,nm1
      am=1.d10
c <<<
      jk=1
c >>>
      do 333 j=1,n
      if(idat(j).ne.0) go to 333
ccxx      if(am.lt.dx(j)) go to 333
         ddd = dabs(am-dx(j))
         if(ddd.ne.0) ddd=ddd/max(dabs(am), dabs(dx(j)))
      if(am.lt.dx(j) .and. ddd.gt.eps01) go to 333

      am=dx(j)
      jk=j
  333 continue
cxx      if(nov.lt.0.and.dx(jk).lt.0.) nw1=i+1
      idat(jk)=i+1
      idata(i+1)=jk
      if(ikn.lt.ite(jk)) ikn=ite(jk)
  332 continue
cxx      if(nov.lt.0) nov=nw1
cc      write(4,2019) (title(i,i3),i=1,20)
cc      write(4,2020)
      as=0.
c <<<
      i41=idata(2)
      as1=dx(i41)
c >>>
      do 335 k2=1,nm1
      i4=idata(k2+1)
c     if(ity(i4).lt.0) go to 335
      if(dx(i4).gt.2.d9) go to 335
      if(k2.ne.1) as=dx(i4)-as
cc      if(k2.eq.1) as1=dx(i4)
cc      as2=dx(i4)-as1
cc      as3=exp(-1./2.*as2)
cc      if(dx(i4).lt.9.d8)
cc     &write(4,2021) k2,(title(i,i4),i=1,20),ite(i4),dx(i4),as2,as3
cc     &               as3(i4)
cc      if(dx(i4).gt.9.d8)
cc     &write(4,2036) k2,(title(i,i4),i=1,20),ite(i4)
      as=dx(i4)
  335 continue
cc      write(4,2029)
cc      write(4,2031)
cc      do 336 k2=1,nm1
cc      i4=idata(k2+1)
cc      if(ity(i4).lt.0) go to 336
cc      if(k2.ne.1) as=dx(i4)-as
cc      if(k2.eq.1) as1=dx(i4)
cc      as2=dx(i4)-as1
cc      as3=exp(-1./2.*as2)
cc      is3=as3*50+0.5
cc      if(is3.ne.0) write(4,2032) k2,(title(i,i4),i=1,10),(st,i=1,is3)
cc      as=dx(i4)
cc  336 continue
cc      write(4,2033)
cc      write(4,2028)
cc      write(4,2017)
      ikpp=ikp
cc      if(ikpp.eq.0) ikpp=11
cc      if(ikpp.gt.n) ikpp=n
cc      if(ikpp.lt.0) ikpp=n
      do 920 k=2,ikpp
cc      if(k.le.10) write(4,2025) k-1
cc      if(k.gt.10.and.k.le.100) write(4,2026) k-1
cc      if(k.gt.100) write(4,2027) k-1
      k3=idata(k)
cc      jj1=noo(k3)+1
cc      jj2=noo(k3)+ite(k3)
      jj1=1
      jj2=ite(k3)
      jj3=ite(k3)
      if(ity(k3).lt.0) go to 920
cc      ii=0
      do 800 j=1,jj3
cc      ttrr(j)=0
      ttrr(k3,j)=0
  800 continue
      do 810 i=i1,i2
cc      jj=0
      do 810 j=jj1,jj2
cc      jj=jj+1
cc      ttrr(jj)=ttrr(jj)+iab(i,j)
      ttrr(k3,j)=ttrr(k3,j)+iab(k3,i,j)
  810 continue
cc      ii=0
      do 820 i=i1,i2
cc      ii=ii+1
cc      jj=0
      do 830 j=jj1,jj2
cc      jj=jj+1
cc      if(ttrr(jj).eq.0) pa(ii,jj)=0.
cc      if(ttrr(jj).eq.0) go to 830
cc      pa(ii,jj)=(iab(i,j)*100.)/ttrr(jj)
      if(ttrr(k3,j).eq.0) pa(k3,i,j)=0.
      if(ttrr(k3,j).eq.0) go to 830
      pa(k3,i,j)=(iab(k3,i,j)*100.)/ttrr(k3,j)
  830 continue
cc      ptc(ii)=(totalc(ii)*100.)/samp
      ptc(i)=(totalc(i)*100.)/samp
  820 continue
cc      jj=0
      do 840 j=jj1,jj2
cc      jj=jj+1
cc      if(ttrr(jj).eq.0) ptr(jj)=0.
cc      if(ttrr(jj).eq.0) go to 840
cc      ptr(jj)=(ttrr(jj)*100.)/ttrr(jj)
      if(ttrr(k3,j).eq.0) ptr(k3,j)=0.
      if(ttrr(k3,j).eq.0) go to 840
      ptr(k3,j)=(ttrr(k3,j)*100.)/ttrr(k3,j)
  840 continue
      pt=(nsamp/samp)*100.
      ip1=1
cc      ip2=ii
      ip2=i2
cc  845 continue
cc      if(ip2.gt.ip1+7) ip2=ip1+7
cc      ncpr=9+(ip2-ip1+2)*4
cc      write(4,2022) (mx,i=1,ncpr)
cc      write(4,2002)(bl,(title(i,i3),i=1,20),j=1,lq)
cc      write(4,2003) (i,i=ip1,ip2)
cc      write(4,2005)(bl,(title(i,k3),i=1,20),j=1,lq)
cc      write(4,2022) (mx,i=1,ncpr)
cc      jj=0
cc      do 890 j=jj1,jj2
cc      jj=jj+1
cc      write(4,2006) jj ,(iab(i,j),i=ip1,ip2),ttrr(jj)
cc  890 continue
cc      write(4,2022) (mx,i=1,ncpr)
cc      write(4,2010) (totalc(i),i=ip1,ip2),nsamp
cc      write(4,2022) (mx,i=1,ncpr)
cc      write(4,2012)
cc      if(ip2.eq.ii) go to 855
cc      ip1=ip2+1
cc      ip2=ii
cc      go to 845
cc  855 continue
cc      ip1=1
cc      ip2=ii
cc  856 continue
cc      if(ip2.gt.ip1+7) ip2=ip1+7
cc      ncpr1=7+(ip2-ip1+2)*6
cc      write(4,2022) (mx,i=1,ncpr1+9)
cc      write(4,2002) bl,(title(i,i3),i=1,20)
cc      write(4,2013) (i,i=ip1,ip2)
cc      write(4,2005) bl,(title(i,k3),i=1,20)
cc      write(4,2022) (mx,i=1,ncpr1+9)
cc      jj=0
cc      do 910 j=jj1,jj2
cc      jj=jj+1
cc      fm1(1)=fm3(1)
cc      fm1(2)=fm3(2)
cc      fm1(4)=fm2(ip2-ip1+1)
cc      write(4,fm1) jj,(pa(i,jj),i=ip1,ip2),ptr(jj),ttrr(jj)
cc  910 continue
cc      write(4,2022) (mx,i=1,ncpr1+9)
cc      fm1(1)=fm3(3)
cc      fm1(2)=fm3(4)
cc      fm1(4)=fm2(ip2-ip1+1)
cc      write(4,fm1) (ptc(i),i=ip1,ip2),pt,nsamp
cc      write(4,2022) (mx,i=1,ncpr1+9)
cc      if(ip2.eq.ii) go to 905
cc      ip1=ip2+1
cc      ip2=ii
cc      go to 856
cc  905 continue
cc      write(4,2024) (title(i,k3),i=1,20)
cc      do 923 j=1,jj
cc      do 923 j=1,jj2
cc      if(it.ne.0.and.xx(k3).ne.0.0) write(4,2008) j,ab(k3,j),
cc     &  ab(k3,j+1)
cc      if(it.eq.0.or.xx(k3).eq.0.0) write(4,2009) j,iaa(k3,j)
cc  923 continue
cc      write(4,2030)
cc  900 continue
  920 continue
cc      novv=0
cc      do 921 i=1,n
cc  921 if(ity(i).ge.0) novv=novv+1
cc      if(izu.eq.1) call pr3(n11,n22,ite,title,ptc,iab,idata,i3,i1,i2,
cc     1                      nx1,xx,it,novv,ikpp)
cc      if(ipart.eq.1) write(4,2001)
      if(ipart.eq.1) go to 315
c
c     storage allocation for necessary variables
c
      nx=n
cc      n=min(nov,novv)
      ikn=ikn+1
cc      ikkk=5000
cc      ikr=1500
cc      ikh=50
cc      ikf=5000
cxx      ikr=4*nx
      ikh=ite(1)
      do 922 i=2,nx
         if(ikh.lt.ite(i)) ikh=ite(i)
  922 continue
c
cc      i11=1+n
cc      i12=i11+ikkk
cc      i13=i12+n
cc      i14=i13+ikkk
cc      i15=i14+ikkk
cc      i16=i15+2*n
cc      i17=i16+ikh*n
cc      i18=i17+n
cc      i19=i18+n
cc      i20=i19+2*n*ikn
cc      i21=i20+n
cc      i22=i21+2*n
cc      i23=i22+n
cc      i24=i23+n
cc      i25=i24+ikr*10
cc      i26=i25+ikr
cc      i27=i26+ikr
cc      i28=i27+20*n
cc      i29=i28+ikr
cc      i30=i29+n
cc      i31=i30+n*ikn
cc      i32=i31+n*ikn
cc      i33=i32+n*ikn
cc      i34=i33+n*ikn
cc      i35=i34+ikr*10
cc      i36=i35+ikr
cc      i37=i36+ikr
cc      i38=i37+20*n
cc      i39=i38+n*nsamp
cc      i40=i39+n
cc      j11=1+ikkk
cc      j12=j11+ikh
cc      j13=j12+ikr
cc      j14=j13+ikr
cc      j15=j14+n
cc      imax1=imax+i40
cc      jmax1=jmax+j15
cc      nnn=ialim-imax1
cc      jn=jalim-jmax1
cc      if(imax1.gt.ialim.or.jmax1.gt.jalim) go to 10
cc      imax0=imax
cc      jmax0=jmax
cc      imax=imax1
cc      jmax=jmax1
cc      if(imax1.gt.imaxx) imaxx=imax1
cc      if(jmax1.gt.jmaxx) jmaxx=jmax1
cc      call mdap(ite,nsamp,iaa,idata(1),idata(i11),idata(i12),idata(i13),
cc     1          idata(i14),idata(i15),idata(i16),ity,idata(i17),
cc     2          idata(i18),idata(i19),idata(i20),idata(i21),idata(i22),
cc     3          idata(i23),idata(i24),idata(i25),idata(i26),idata(i27),
cc     4          idata(i28),idata(i29),idata(i30),idata(i31),idata(i32),
cc     5          idata(i33),idata(i34),idata(i35),idata(i36),idata(i37),
cc     6          idata(i38),idata(i39),tttr(j1),tttr(j11),tttr(j12),
cc     7          tttr(j13),tttr(j14),ab,xx,iconv,title,item1,idata(i40),
cc     8          tttr(j15),it,ikkk,ikr,ikh,ikn,ikf,n,n33,recode,nx,nnn,
cc     9          jn,izu,icl,icls)
      call mdap0(ite,nsamp,iaa,idata,ida,iby,ibc,aic1,caa,icaa,nnia,pbc,
     1           aaam,da,ab,iabse,xx,iconv,item1,ity,it,ikkk,ikr,jkr,
     2           ikh,ikn,lk77,morder,nov,n11,n33,recode,nx,
     3           icl,icls,eps01,ier)

      n = nx

  315 continue
cc      if(ik.eq.l) go to 310
cc      i3=face(ik+1)
cc      i2=i2+item(i3)
  310 continue
c>>>
      n = nx1
      if(ier(1).eq.2048) ier(3) = ikh
c>>>
      return
cc   10 write(4,2016) imax1,ialim,jmax1,jalim
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
cc  930 write(4,2034) (title(i,i3),i=1,20)
  930 ier(1)=2034
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
      return
cc  940 write(4,2035) i,k
  940 ier(1)=2035
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
      return
 2001 format(1h )
 2002 format(8x,a1,'(',20a1,')',37x,a1,'(',20a1,')')
 2003 format(7x,20i5)
 2004 format(1h+,60x,7x,10i6)
 2005 format(a1,'(',20a1,')',37x,a1,'(',20a1,')')
 2006 format(i5,2x,20i5)
 2007 format(1h+,60x,i5,2x,10f6.1)
 2008 format(i10,15x,d12.5,' - ',d12.5)
 2009 format(i10,20x,i10)
 2010 format('total',2x,20i5)
 2011 format(1h+,60x,'total',2x,10f6.1)
 2012 format(' ')
 2013 format(7x,20i6)
 2014 format(i5,2x,20f6.1)
 2015 format('total',2x,20f6.1)
 2016 format(' ia or a dimension over ',4i10)
 2017 format('+',41('-'),'+'/
     1       'i the coresponding two-way tables         i'/
     2       'i (each having an optimal categorization) i'/
     3       '+',41('-'),'+')
 2018 format(' response variable  : ',20a1)
 2019 format(72('=')//'<2>'/'+',38('-'),'+'/
     1       'i list of single explanatory variables i'/
     2       'i (arranged in ascending order of aic) i'/'+',38('-'),'+'/
     3       'response variable  : ','(',20a1,')')
 2020 format(72('-')/
     1       '  no.',2x,'explanatory',    4x,'number of categories',
     2 3x,'    ',4x,'difference       '/
     3 7x,'variable  ',5x,'of exp. var. ',7x,'  a i c ',2x,' of aic',
     4 6x,'weight'/72('-'))
 2021 format(i5,3x,20a1,1x,i5,5x,f10.2,1x,f10.2,1x,f10.5)
 2022 format(100a1)
 2023 format(1h+,60x,60a1)
 2024 format('<note>'/5x,20a1,5x,'class interval')
 2025 format('<3',i1,'>')
 2026 format('<3',i2,'>')
 2027 format('<3',i3,'>')
 2028 format(72('=')//'<3>')
 2029 format(72('-'))
 2030 format(72('+'))
 2031 format( ' weight'/
     1       15x,'0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0'
     2       /16x,10('+----'),'+')
 2032 format(i4,1x,10a1,1x,'i',50a1)
 2033 format(16x,10('+----'),'+')
 2034 format(//' error'/
     1       ' the response variable "',20a1,'" is continuous.'/
     2       ' a continuous variable cannot be specified as the res',
     3       'ponse variable.')
 2035 format(//' error'/
     1       ' the value of ',i3,'-th variable of ',i5,'-th ',
     2       'unit'/' is beyond the interval specified in cols 26-28',
     3       ' and 29-30.')
 2036 format(i5,3x,20a1,1x,i5,'(1)',6x,3(' - - -   ',1x))
      end
cc      subroutine mdap(item,nsamp,iaa,imm,ibd,iby,ib,ibc,items,itemz,ity,
cc     1                idt,itemx,ias1,itemy,itemt,lc,lca,ca,ica,lg,cb,
cc     2                nni,itype,ias2,icon,icon1,icon2,caa,icaa,nnia,c,
cc     3                data,idata,pbc,aicmm,aa,aaa,da,ab,xx,iconv,title,
cc     4                item1,iw,w,itx,ikkk,ikr,ikh,ikn,ikf,n,n33,recode,
cc     5                nx,nnn,jnn,izu,icl,icls,in)
      subroutine mdap0(item,nsamp,iaa,imm,idata,iby,ibc,aic1,caa,icaa,
     1              nnia,pbc,aaa,da,ab,iabse,xx,iconv,item1,ity,itx,
     2              ikkk,ikr,jkr,ikh,ikn,lk77,morder,n,n11,n33,recode,
     3              nx,icl,icls,eps01,ier)
c
c     this subroutine searches for the optimal multidimensional
c     contingency table with respect to the combination of explanatory
c     variables and the categorization of each one.
c
      implicit real*8(a-h,o-z)
cc      integer*2 ib,ias1,itemy,lc,itemt,c,lca,ca,ica,lg,cb,ibc,items,
cc     1          ity,itype,ias2,iw,iaa,idt,nni,data,imm,icon,icon1,
cc     2          icon2,caa,icaa,nnia,idata,iconv,title,item1,itemz,
cc     3          item,itemx,ibd,iby,icls
cc      integer*2 title
      integer ca,data,caa
      integer recode
cc      dimension ibd(ikkk),iby(n),ib(ikkk),pbc(ikkk),ibc(ikkk),
cc     1          items(2,n),itemz(ikh,n),c(20,n),idt(n),idata(nx),
cc     2          item(nx),itemx(n),ias1(2,n,ikn),ity(nx),icon2(n,ikn),
cc     3          aicmm(ikh),itemy(n),itemt(n,2),lc(n),lca(n),da(nx),
cc     4          aa(ikr),ca(ikr,10),ica(ikr),lg(ikr),cb(20,n),nni(ikr),
cc     5          itype(n),ias2(n,ikn),ab(nx,n33),iaa(nx,n33),iw(nnn),
cc     6          w(jnn),data(nsamp,n),imm(n),icon(n,ikn),icon1(n,ikn),
cc     7            ica(ikr),lg(ikr),nni(ikr),itype(n),ias2(n,ikn),
cc     8            data(nsamp,n),icon(n,ikn),icon1(n,ikn),
cc     9            icaa(ikr),nnia(ikr),pbcc(n11),morder(ikr)
      dimension iaa(nx,n33),imm(nx),idata(nx,nsamp),aic1(icl+1),
     1           iby(n,ikkk,icl+1),ibc(n11,ikkk,icl+1),
     2           pbc(n11,ikkk,icl+1),caa(ikr,jkr),icaa(ikr),
     3           nnia(ikr),aaa(ikr),da(nx,nsamp),ab(nx,n33),xx(nx),
     4           iconv(20,recode),item1(nx),morder(ikr),
     4           icls(max(10,n+1),icl),
     5           ibd(n11),ib(ikkk),item(nx),itemx(n),items(2,n),ity(nx),
     6           itemt(n,2),itemy(n),itemz(ikh,n),ias1(2,n,ikn),
     7           ias2(n,ikn),lca(n),aa(ikr),aicmm(ikh),pbcc(n11),lc(n),
     8           ca(ikr,jkr),ica(ikr),lg(ikr),nni(ikr),itype(n),idt(n),
     9           data(nsamp,n),icon(n,ikn),icon1(n,ikn),icon2(n,ikn),
     &           ier(2),iabse(n,n33,icl+1)
cc      character*4 fmt(20),fm1(5),fm2(5),fm3(4),fm4(10),fm11(10)
cc      character*4 fm12(10),fm13(10),fm21(10),fm22(10),fm23(10)
cc      character*1 mxx(6)
cc      character*2 cz(20)
cc      common ialim,imax,jalim,jmax,imaxx,jmaxx
cc      data mxx/'x','(',')',':','-',' '/
cc      data (cz(i),i=1,20)
cc     &       /'  ','- ','  ','- ','  ','- ',14*'  '/
cc      data (fmt(i),i=1,20)
cc     &     /'(1x,','    ','    ','i1, ','    ','20i5',')   ',13*'    '/
cc      data (fm1(i),i=1,5)
cc     &        /'   1','   2','   3','   4','   5'/
cc      data (fm2(i),i=1,5)
cc     &        /'14x,','13x,','12x,','11x,','10x,'/
cc      data (fm3(i),i=1,4)
cc     &        /'10i5',')   ','10f6','.1) '/
cc      data (fm4(i),i=1,10)
cc     &        /'f6.1',',2h ','(,i5',',1h)',' )  ',5*'    '/
cc      data (fm11(i),i=1,10)
cc     &      /'(1x,','    ','i1, ','    ','5x, ','10i5',')   ',3*'    '/
cc      data (fm12(i),i=1,10)
cc     &  /'   1','   2','   3','   4','   5','   6','   7','   8','   9',
cc     &   '    '/
cc      data (fm13(i),i=1,10)
cc     &  /'9x, ','8x, ','7x, ','6x, ','5x, ','4x, ','3x, ','2x, ','1x, ',
cc     &   '    '/
cc      data (fm21(i),i=1,10)
cc     &         /'(1x,','    ','a1, ','    ','5x, ','55a1',')   ',
cc     &         3*'    '/
cc      data (fm22(i),i=1,10)
cc     &  /'   1','   2','   3','   4','   5','   6','   7','   8','   9',
cc     &   '    '/
cc      data (fm23(i),i=1,10)
cc     &  /'9x, ','8x, ','7x, ','6x, ','5x, ','4x, ','3x, ','2x, ','1x, ',
cc     &   '    '/
c
c     discretization of original data according to the optimal
c     categorization of each explanatory variable and store
c     of the data
c
cc      icl=icl-1
cc      rewind 1
      i3=imm(1)
      if(itx.ne.0)    go to 920
      kk5=item(i3)
      do 930 i=1,kk5
      iaa(i3,i)=i+item1(i3)-1
  930 continue
  920 continue
cc      do 940 j=1,20
cc      c(j,1)=title(j,i3)
cc  940 continue
      ii=1
cc      do 950 i=2,n
cc      i3=imm(i)
cc      do 960 j=1,20
cc      c(j,i)=title(j,i3)
cc  960 continue
cc  950 continue
      do 990 i=1,nsamp
cc      if(itx.eq.0) read(1) (idata(j),j=1,nx)
cc      if(itx.ne.0) read(1)(da(j),j=1,nx)
      do 1100 iy=1,n
      iy3=imm(iy)
c <<<
      j=1
c >>>
cc      if(itx.ne.0.and.xx(iy3).eq.0) j=da(iy3)+0.1
      if(itx.ne.0.and.xx(iy3).eq.0) j=da(iy3,i)+0.1
      if(itx.ne.0.and.xx(iy3).eq.0) go to 1110
      kk5=item(iy3)
      do 1120 j=1,kk5
      if(itx.ne.0.and.xx(iy3).ne.0.) go to 1130
cc      if(idata( iy3).le.iaa(iy3,j)-item1(iy3)+1) go to 1110
      if(idata( iy3,i).le.iaa(iy3,j)-item1(iy3)+1) go to 1110
      go to 1120
cc 1130 if(abs(da(iy3)-ab(iy3,j+1)).lt.xx(iy3)/1.d5) go to 1120
cc      if(ab(iy3,j+1).gt.da(iy3)) go to 1110
 1130 if(abs(da(iy3,i)-ab(iy3,j+1)).lt.xx(iy3)/1.d5) go to 1120
      if(ab(iy3,j+1).gt.da(iy3,i)) go to 1110
 1120 continue
 1110 continue
      if(j.gt.item(iy3).and.xx(iy3).ne.0.) j=item(iy3)
      data(i,iy)=j
 1100 continue
      if(recode.eq.0.or.itx.eq.0) go to 1240
      do 1250 iy=1,n
      ig1=imm(iy)
      do 1260 iz=1,recode
      if(iconv(1,iz).ne.ig1) go to 1260
      ik1=data(i,iy)+1
      if(data(i,iy).eq.0) ik1=11
      data(i,iy)=iconv(ik1,iz)
      go to 1250
 1260 continue
 1250 continue
 1240 continue
      if(itx.eq.0) go to 1280
      do 1270 iy=1,n
      iy3=imm(iy)
      if(xx(iy3).ne.0.) go to 1270
      kk5=item(iy3)
      do 1290 j=1,kk5
      if(data(i,iy).le.iaa(iy3,j)) go to 1300
 1290 continue
      j=kk5
 1300 data(i,iy)=j
 1270 continue
 1280 continue
      do 1275 iy=1,n
      iy3=imm(iy)
      kk5=item(iy3)
      if(data(i,iy).gt.kk5) data(i,iy)=kk5
 1275 continue
  990 continue
      do 10 i=1,n
      ii=imm(i)
      itemx(i)=item(ii)
      itemy(i)=item(ii)
      items(1,i)=item(ii)
      itemz(1,i)=item(ii)
      itype(i)=ity(ii)
   10 continue
      do 15 i=1,ikr
      aa(i)=1.d10
   15 continue
      do 16 i=1,ikh
      aicmm(i)=1.d10
   16 continue
      do 20 i=1,n
      itemm=itemx(i)+1
      do 20 j=1,itemm
      ias1(1,i,j)=1
      ias2(i,j)=1
      if(j.ne.itemm) icon(i,j)=j
      if(j.ne.itemm) icon1(i,j)=j
   20 continue
c
c     storage allocation for the necessary variables
c
cc      i10=1+n
cc      i11=i10+n*n
cc      i12=i11+n*n
cc      i13=i12+n
cc      i14=i13+n
cc      i15=i14+n
cc      i16=i15+n
cc      i17=i16+2*n
cc      i18=i17+2
cc      i19=i18+10*n
cc      i20=i19+n
cc      j11=1+n
cc      j12=j11+n
cc      j13=j12+n
cc      j14=j13+n
cc      imax1=imax+i20
cc      jmax1=jmax+j14
cc      nnn=ialim-imax1
cc      if(imax1.gt.ialim.or.jmax1.gt.jalim) go to 570
cc      imax0=imax
cc      jmax0=jmax
cc      imax=imax1
cc      jmax=jmax1
cc      if(imax1.gt.imaxx) imaxx=imax1
cc      if(jmax1.gt.jmaxx) jmaxx=jmax1
c
c     stepwise selection procedure of explanatory variables
c
cc      call aicp(n,itemx,lc,lk2,aic,iw(i10),iw(i11),iw(i12),iw(i13),
cc     1          iw(i14),iw(i15),iw(i16),iw(i17),c,iw(i18),iw(i19),
cc     2          ca,ica,nni,data,icon,w(1),w(j11),w(j12),w(j13),aa,
cc     3          lk7,ikr,ikn,ikf,n,nsamp,iw(i20),nnn,in)
      call aicp0(n,itemx,lc,lk2,aic,ca,ica,nni,data,icon,aa,lk7,ikr,jkr,
     1          ikn,n,nsamp,eps01,ier)
      if( ier(1).ne.0 ) return

      aicmm(1)=aic
      kt=1
      acmmm=aic
      lk5=lk2
      lk77=lk7
      do 40 i=1,lk5
   40 lca(i)=lc(i)
      do 35 i=1,lk77
      aaa(i)=aa(i)
      icaa(i)=ica(i)
      nnia(i)=nni(i)
         if(icaa(i).gt.jkr) go to 573
cc      do 35 ii=1,10
      do 35 ii=1,icaa(i)
      caa(i,ii)=ca(i,ii)
   35 continue
      do 31 i=1,n
      itemxx=itemz(1,i)
      do 30 iii=1,itemxx
      icon2(i,iii)=icon(i,iii)
   30 continue
   31 continue
c
c     pooling of the categories of each explanatory variable in a
c     multidimensional table
c
   60 continue
      kt=kt+1
      aicm=1.d10
      jj=0
c <<<
      ll=1
c >>>
      do 65 i=1,n
      do 65 iii=1,ikn
      icon(i,iii)=icon1(i,iii)
   65 continue
      do 70 j=2,n
      if(itemx(j).le.2) go to 70
      if(itype(j).eq.1) go to 71
      if(itype(j).eq.2) go to 70
c
c     pooling of categories of each explanatory variable in a
c     multidimensional table in the case of 'ity(j)=0'
c
      do 106 i=1,n
      do 106 l=1,2
      itemt(i,l)=itemx(i)
      if(l.eq.1.and.i.eq.j) itemt(i,l)=(itemx(i)+1)/2
      if(l.eq.2.and.i.eq.j) itemt(i,l)=itemx(i)/2+1
  106 continue
      do 130 l=1,2
      do 133 i=1,n
      itemxx=itemz(1,i)
      do 133 iii=1,itemxx
      icon(i,iii)=icon1(i,iii)
  133 continue
      itemxx=itemz(1,j)
      do 135 i=1,itemxx
      if(l.eq.1) icon(j,i)=(icon1(j,i)+1)/2
      if(l.eq.2) icon(j,i)=icon1(j,i)/2+1
  135 continue
cc      imax=imax1
cc      jmax=jmax1
cc      call aicp(n,itemt(1,l),lc,lk2,aic,iw(i10),iw(i11),iw(i12),
cc     1          iw(i13),iw(i14),iw(i15),iw(i16),iw(i17),c,iw(i18),
cc     2          iw(i19),ca,ica,nni,data,icon,w(1),w(j11),w(j12),
cc     3          w(j13),aa,lk7,ikr,ikn,ikf,n,nsamp,iw(i20),nnn,in)
      call aicp0(n,itemt(1,l),lc,lk2,aic,ca,ica,nni,data,icon,aa,lk7,
     1          ikr,jkr,ikn,n,nsamp,eps01,ier)
      if( ier(1).ne.0 ) return

ccxx      if(aicm.lt.aic) go to 130
         daic=dabs(aicm-aic)
         if(daic.ne.0) daic=daic/max(dabs(aicm), dabs(aic))
      if(aicm.lt.aic .and. daic.gt.eps01) go to 130

      aicm=aic
      iko=ik2
      jj=j
      ll=l
      do 150 i=1,n
      itemy(i)=itemx(i)
      if(l.eq.1.and.i.eq.j) itemy(j)=(itemx(j)+1)/2
      if(l.eq.2.and.i.eq.j) itemy(j)=itemx(j)/2+1
  150 continue

ccxx      if(acmmm.lt.aic+1.d-10) go to 160
      if(acmmm.le.aic) go to 160
        daic=dabs(acmmm-aic)
        if(daic.ne.0) daic=daic/max(dabs(acmmm), dabs(aic))
        if(daic.lt.eps01) go to 160

      acmmm=aic
      lk5=lk2
      lk77=lk7
      do 170 i=1,lk5
  170 lca(i)=lc(i)
      do 155 i=1,lk77
      aaa(i)=aa(i)
      icaa(i)=ica(i)
      nnia(i)=nni(i)
         if(icaa(i).gt.jkr) go to 573
cc      do 155 ii=1,10
      do 155 ii=1,icaa(i)
      caa(i,ii)=ca(i,ii)
  155 continue
c <<<
      do 168 i=1,n
      itemx1=itemz(1,i)
      do 162 iii=1,itemx1
      icon2(i,iii)=icon(i,iii)
  162 continue
  168 continue
c >>>
  160 continue
  130 continue
      go to 70
   71 continue
c
c     the similar processing in the case of  'ity(j)=1'
c
      itemxx=itemx(j)-1
      do 73 ii=1,itemxx
      itemx(j)=itemxx
      iii=0
      itemx1=itemz(1,j)
      do 74 i=1,itemx1
      if(ii+1.gt.icon1(j,i)) icon(j,i)=icon1(j,i)
      if(ii+1.le.icon1(j,i)) icon(j,i)=icon1(j,i)-1
   74 continue
cc      imax=imax1
cc      jmax=jmax1
cc      call aicp(n,itemx,lc,lk2,aic,iw(i10),iw(i11),iw(i12),iw(i13),
cc     1          iw(i14),iw(i15),iw(i16),iw(i17),c,iw(i18),iw(i19),
cc     2          ca,ica,nni,data,icon,w(1),w(j11),w(j12),w(j13),aa,
cc     3          lk7,ikr,ikn,ikf,n,nsamp,iw(i20),nnn,in)
      call aicp0(n,itemx,lc,lk2,aic,ca,ica,nni,data,icon,aa,lk7,ikr,jkr,
     1          ikn,n,nsamp,eps01,ier)
      if( ier(1).ne.0 ) return

      itemx(j)=itemx(j)+1
      if(aicm.le.aic) go to 73
         daic=dabs(aicm-aic)
         if(daic.ne.0) daic=daic/max(dabs(aicm), dabs(aic))
         if(daic.lt.eps01) go to 73

      aicm=aic
      iko=ik2
      jj=j
      ll=ii
      itemy(j)=itemxx
      if(acmmm.le.aic) go to 85
         daic=dabs(acmmm-aic)
         if(daic.ne.0) daic=daic/max(dabs(acmmm), dabs(aic))
         if(daic.lt.eps01) go to 85

      acmmm=aic
      lk5=lk2
      lk77=lk7
      do 81 i=1,lk5
   81 lca(i)=lc(i)
      do 185 i=1,lk77
      aaa(i)=aa(i)
      icaa(i)=ica(i)
      nnia(i)=nni(i)
         if(icaa(i).gt.jkr) go to 573
cc      do 185 iii=1,10
      do 185 iii=1,icaa(i)
      caa(i,iii)=ca(i,iii)
  185 continue
      do 88 i=1,n
      itemx1=itemz(1,i)
      do 82 iii=1,itemx1
      icon2(i,iii)=icon(i,iii)
   82 continue
   88 continue
   85 continue
   73 continue
      do 83 i=1,ikn
   83 icon(j,i)=icon1(j,i)
   70 continue
c
c     keeping the categorization with a temporary maice
c
      if(jj.eq.0) go to 275
      aicmm(kt)=aicm
      ikk=iko
      items(2,jj)=items(1,jj)
      items(1,jj)=itemy(jj)
      do 200 i=1,n
      itemz(kt,i)=items(1,i)
  200 continue
      itemmx=itemx(jj)+1
      do 210 i=1,itemmx
      ias1(2,jj,i)=ias1(1,jj,i)
  210 continue
      ias1(1,jj,1)=ias1(2,jj,1)+1
      itemmy=itemy(jj)
      itemxx=itemz(1,jj)
      if(itype(jj).eq.1) go to 236
      do 235 i=1,itemxx
      if(ll.eq.1) icon1(jj,i)=(icon1(jj,i)+1)/2
      if(ll.eq.2) icon1(jj,i)=icon1(jj,i)/2+1
  235 continue
      go to 237
  236 continue
      do 238 i=1,itemxx
      if(ll+1.gt.icon(jj,i)) icon1(jj,i)=icon(jj,i)
      if(ll+1.le.icon(jj,i)) icon1(jj,i)=icon(jj,i)-1
  238 continue
  237 continue
      if(itype(jj).eq.1) go to 245
      if(ll.eq.2) go to 220
      ii=0
      do 230 i=1,itemmy
      if(itemx(jj).lt.ii+1) ias1(2,jj,ii+2)=0
      ias1(1,jj,i+1)=ias1(2,jj,ii+1)+ias1(2,jj,ii+2)
      ii=ii+2
  230 continue
      go to 240
  220 continue
      ias1(1,jj,2)=ias1(2,jj,2)
      ii=2
      do 250 i=2,itemmy
      if(itemx(jj).lt.ii  ) ias1(2,jj,ii+1)=0
      ias1(1,jj,i+1)=ias1(2,jj,ii)+ias1(2,jj,ii+1)
  250 continue
  240 continue
      go to 246
  245 continue
      do 256 i=1,n
      itemy(i)=itemx(i)
      if(jj.eq.i) itemy(i)=itemy(i)-1
  256 continue
      ii=1
      ias1(1,jj,ii)=ias1(2,jj,1)
      do 255 i=2,itemmy
      if(ll.ne.i+1) ii=ii+1
      ias1(1,jj,ii)=ias1(2,jj,i)
      if(ll.ne.i+1) ias1(1,jj,ii)=ias1(1,jj,ii)+ias1(2,jj,i)
  255 continue
  246 continue
ccxx      if(acmmm.lt.aicm) goto 252
         daic=dabs(acmmm-aicm)
         if(daic.ne.0) daic=daic/max(dabs(acmmm), dabs(aicm))
      if(acmmm.lt.aicm .and. daic.gt.eps01) goto 252

      itemmy=itemy(jj)
      do 251 i=1,itemmy
      ias2(jj,i)=ias1(1,jj,i+1)
  251 continue
  252 continue
      do 260 i=1,n
      itemx(i)=itemy(i)
  260 continue
      if(itype(jj).eq.1) go to 270
      if(kt.le.3) go to 270
cc      imax=imax0
cc      jmax=jmax0
cc      i11=1+n
cc      imax=imax+n
c
c     re-dividing the pooled categories
c
cc      call bun(icon,ias1,jj,aic,items,ix,itemy,lk2,iw(1),ca,ica,nni,c,
cc     1         aa,n,iw(i11),w(1),data,nsamp,ikr,ikn,ikf,lk7,nnn,in)
      call bun0(icon,ias1,jj,aic,items,ix,itemy,lk2,ca,ica,nni,aa,n,
     1         data,nsamp,ikr,jkr,ikn,lk7,eps01,ier)
      if(ier(1).ne.0) return

      if(aicmm(kt-1).le.aic.or.aicmm(kt).le.aic) go to 270
ccxx      if(aic-aicmm(kt-1).gt.1.d-10.or.aic-aicmm(kt).gt.1.d-10)
      if(aic-aicmm(kt-1).gt.epa01 .or. aic-aicmm(kt).gt.eps01)
     1 go to 270
      kt=kt+1
      aicmm(kt)=aic
      do 280 i=1,n
      itemx(i)=itemy(i)
  280 continue
      itemxx=itemz(1,ix)
      do 335 i=1,itemxx
      icon1(ix,i)=icon(ix,i)
  335 continue
ccxx      if(acmmm.lt.aic) go to 290
         daic=dabs(acmmm-aic)
         if(daic.ne.0) daic=daic/max(dabs(acmmm), dabs(aic))
      if(acmmm.lt.aic .and. daic.gt.eps01) go to 290

      acmmm=aic
      lk5=lk2
      lk77=lk7
cc      do 305 i=1,200
      do 305 i=1,lk77
      aaa(i)=aa(i)
      icaa(i)=ica(i)
      nnia(i)=nni(i)
         if(icaa(i).gt.jkr) go to 573
cc      do 305 ii=1,10
      do 305 ii=1,icaa(i)
      caa(i,ii)=ca(i,ii)
  305 continue
      do 300 i=1,n
      itemxx=itemz(1,i)
      do 300 iii=1,itemxx
      icon2(i,iii)=icon(i,iii)
  300 continue
      do 320 i=1,n
      itemz(kt,i)=itemx(i)
  320 continue
  290 continue
      items(1,ix)=items(2,ix)
      itemyy=items(2,ix)
      do 330 i=1,itemyy
      ias1(1,ix,i)=ias1(2,ix,i)
  330 continue
      do 331 i=1,itemyy
      ias2(ix,i)=ias1(1,ix,i+1)
  331 continue
  270 continue
      if(aicmm(kt).lt.aicmm(kt-1).and.aicmm(kt-1)-aicmm(kt).gt.1.d-10)
     1  go to 60
c
c     printing out 'aic's at each dimension ' and 'summary of the
c     subsets'
c
  275 continue
      lk=lk77
cc      write(4,2028)
      ktt=2
      do 346 i=1,lk
      if(ktt.lt.icaa(i)) ktt=icaa(i)
  346 continue
      ktt=ktt+1
      do 348 iii=2,ktt
cc      if(iii.ne.ktt) write(4,2035) iii-1
cc      if(iii.eq.ktt) write(4,2036)
cc      if(iii.eq.ktt) write(4,2029)
cc      write(4,2027) (c(ix,1),ix=1,20)
      iiii=iii-1
cc      if(iii.ne.ktt) write(4,2024) iiii
      lkt=0
      do 349 i=1,lk
      if(icaa(i).eq.iii.or.iii.eq.ktt) lkt=lkt+1
  349 continue
cc      write(4,2026)
      do 341 i=1,lk
      lg(i)=0
  341 continue
      aaaa=-1.d10
      iij=0
c <<<
      lk44=2
c >>>
      do 342 ij=1,lk
      aminx=10.**10
cxx      do 343 i=1,lk
      if (lk.lt.2) go to 3430
      do 343 i =2,lk
      if(lg(i).ne.0) go to 343
      if(aminx.le.aaa(i)) go to 343
         daic=dabs(aminx-aaa(i))
         if(daic.ne.0) daic=daic/max(dabs(aminx), dabs(aaa(i)))
         if(daic.lt.eps01) go to 343

      aminx=aaa(i)
      lk3=i
  343 continue
 3430    continue
      lg(lk3)=lk3
      lk4=icaa(lk3)
      if(lk4.ne.iii.and.iii.ne.ktt) go to 342
cc      do 347 i=1,lk4
cc      lx=caa(lk3,i)
cc      do 347 ix=1,20
cc      cb(ix,i)=c(ix,lx)
cc  347 continue
ccxx      if(aaaa.lt.aaa(lk3)-1.0d-10) go to 345
         daic=dabs(aaaa-aaa(lk3))
         if(daic.ne.0) daic=daic/max(dabs(aaaa), dabs(aaa(lk3)))
       if(aaaa.lt.aaa(lk3) .and. daic.gt.eps01) go to 345

cxx      if(lk4.ne.lk44) go to 345
      if(iii.ne.ktt .and. lk4.ne.lk44) go to 345

      if(lk4.lt.2) go to 342
         if(lk4.gt.jkr) go to 574
      call eqck(caa,ikr,lk4,lk3,lk33,ijk)
      if(ijk.eq.0) go to 342
  345 continue
      aaaa=aaa(lk3)
cc      if(iij.eq.0) aaa1=aaa(lk3)
cc      aaa2=aaa(lk3)-aaa1
cc      aaa3=exp(-1./2.*aaa2)
      iij=iij+1
      if(iij.eq.1.and.iii.eq.ktt.and.lk4.eq.1) lk5=0
      if(iij.gt.100) go to 344
c <<<
      if( iii.ne. ktt ) morder(ij)=lk3
c >>>
cc      if(lk4.ge.2) write(4,2018) iij,(cb(ix,2),ix=1,20),nnia(lk3),
cc     1                           aaa(lk3),aaa2,aaa3
cc      if(lk4.ge.3) write(4,2019) ((cb(ix,i),ix=1,20),i=3,lk4)
cc      if(lk4.lt.2) write(4,2018) iij,cz,nnia(lk3),aaa(lk3),aaa2,aaa3
      lk33=lk3
      lk44=lk4
  342 continue
  344 continue
cc      write(4,2042)
  348 continue
c
c     printing out 'contingency table with the optimal combination
c     and categorization of explanatory variables'
c
cc      imax0=imax
cc      jmax0=jmax
      aicm=1.d10
c <<<
      ji=1
c >>>
      do 350 i=1,kt
ccxx      if(aicmm(i).ge.aicm) go to 350
      if(aicm.le.aicmm(i)) go to 350
         daic=dabs(aicm-aicmm(i))
         if(daic.ne.0) daic=daic/max(dabs(aicm), dabs(aicmm(i)))
         if(daic.lt.eps01) go to 350

      aicm=aicmm(i)
      ji=i
  350 continue
      ncount=1
      if(lk5.eq.0) go to 361
  600 continue
cc      imax=imax0
cc      jmax=jmax0
      do 360 i=1,lk5
      lc(i)=lca(i)
      ix=lca(i)
cc      do 365 j=1,20
cc      cb(j,i)=c(j,ix)
cc  365 continue
  360 continue
  361 continue
      lk6=lk5+1
cc      if(ncount.eq.1) write(4,2037)
cc      if(ncount.eq.1) write(4,2005)
cc      if(ncount.ne.1) write(4,2047)
      i=1
cc      write(4,2006) mxx(1),mxx(2),i,mxx(3),mxx(4),(c(j,1),j=1,20),
cc     1                           (mxx(1),mxx(2),i,mxx(3),mxx(4),
cc     1                           (cb(j,i-1),j=1,20),i=2,lk6)
cc      write(4,2002)
cc      if(lk5.ne.0) write(4,2008) (mxx(1),j=2,lk6),(mxx(6),j=1,11-lk6)
      idt(1)=itemz(ji,1)
      ikk=idt(1)
      if(lk5.eq.0) go to 375
      do 370 iii=1,lk5
      lyy=lc(iii)
      idt(iii+1)=itemz(ji,lyy)
      ikk=ikk*idt(iii+1)
  370 continue
      if(ikkk.lt.ikk) go to 575
  375 continue
      do 380 i=1,n
      itemx(i)=itemz(ji,i)
  380 continue
      do 390 i=1,ikk
      ib(i)=0
  390 continue
      do 400 i=1,nsamp
      iy1=1
      iyy=data(i,iy1)
      ite=idt(1)
      iax=icon2(iy1,iyy)
      ik2=ite
      if(lk5.eq.0) go to 425
      do 420 k=2,lk6
      nk1=lk6-k+1
      iy1=lc(nk1)
      ite=idt(nk1+1)
      iyy=data(i,iy1)
      iyy=icon2(iy1,iyy)
      iax=iax+(iyy-1)*ik2
      ik2=ik2*ite
  420 continue
  425 continue
      ib(iax)=ib(iax)+1
  400 continue
      kk=1
      do 430 i=1,lk6
      kk=kk*idt(i)
  430 continue
      do 440 ill=1,2
      kkj=kk/idt(1)
      kx1=idt(1)
      ip1=1
      ip2=kx1
  445 continue
      if(ip2.gt.ip1+7) ip2=ip1+7
      do 450 it=1,kx1
      ibd(it)=0
  450 continue
cc      if(ill.eq.1) ncpr=2+5*(ip2-ip1+2)
cc      if(ill.eq.2) ncpr=2+6*(ip2-ip1+2)
cc      if(lk6.eq.1) go to 455
cc      fm21(2)=fm22(lk6-1)
cc      fm21(4)=fm23(lk6-1)
cc      if(ill.eq.1) write(4,fm21) (mxx(5),i=2,lk6),(mxx(5),i=1,ncpr)
cc      if(ill.eq.2) write(4,2001) (mxx(5),i=1,ncpr)
cc      fm11(2)=fm12(lk6-1)
cc      fm11(4)=fm13(lk6-1)
cc      if(ill.eq.1) write(4,fm11) (it,it=2,lk6),(it,it=ip1,ip2)
cc      if(ill.eq.2) write(4,2033) (it,it=ip1,ip2)
cc      if(ill.eq.1) write(4,fm21) (mxx(5),i=2,lk6),(mxx(5),i=1,ncpr)
cc      if(ill.eq.2) write(4,2001) (mxx(5),i=1,ncpr)
cxx  455 do 460 i=1,kkj
  455 continue
cxx        if(kkj. gt. n33*n33*n33) go to 585
         if(kkj. gt. ikkk) go to 585
      do 460 i=1,kkj
      i01=i-1
      kx=idt(lk6)
cc      iby(lk6)=mod(i01,kx)+1
      iby(lk6,i,ncount)=mod(i01,kx)+1
      ky=kx
      if(lk5.eq.0) go to 475
      do 470 j=2,lk6
      i02=i01/ky
      lk6j=lk6+1-j
      kx=idt(lk6j)
cc      iby(lk6j)=mod(i02,kx)+1
      iby(lk6j,i,ncount)=mod(i02,kx)+1
      ky=ky*kx
  470 continue
  475 continue
      ibct=0
      iq=(i-1)*kx1
      do 480 it=1,kx1
      iq=iq+1
      ibct=ib(iq)+ibct
cc      ibc(it,i)=ib(iq)
      ibc(it,i,ncount)=ib(iq)
      ibd(it)=ibd(it)+ib(iq)
  480 continue
      if(ill.eq.1) go to 500
      do 510 it=1,kx1
cc      pbc(it)=0.
      pbc(it,i,ncount)=0.
      if(ibct.eq.0) go to 510
cc    pbc(it)=ibc(it)*100./ibct
      pbc(it,i,ncount)=ibc(it,i,ncount)*100./ibct      
  510 continue
cc      if(ibct.eq.0) pbct=0.
cc      if(ibct.ne.0) pbct=ibct*100./ibct
  500 continue
cc      if(lk5.eq.0) go to 460
cc      fmt(3)=fm1(lk6-1)
cc      fmt(5)=fm2(lk6-1)
cc      fmt(6)=fm3(2*ill-1)
cc      fmt(7)=fm3(2*ill)
cc      if(ill.eq.2) fmt(6)=fm22(ip2-ip1+2)
cc      do 505 j=1,5
cc      if(ill.eq.2) fmt(j+6)=fm4(j)
cc  505 continue
cc      if(ill.eq.1) write(4,fmt)(iby(j),j=2,lk6),(ibc(it),it=ip1,ip2),
cc     1                         ibct
cc      if(ill.eq.2) write(4,fmt)(iby(j),j=2,lk6),(pbc(it),it=ip1,ip2),
cc     1                          pbct,ibct
  460 continue
      if(ill.eq.1) go to 520
      do 530 it=1,kx1
cc      pbc(it)=ibd(it)*100./nsamp
      pbcc(it)=ibd(it)*100./nsamp      
  530 continue
  520 continue
cc      write(4,2001) (mxx(5),i=1,ncpr)
cc      if(ill.eq.1) write(4,2015) (ibd(it),it=ip1,ip2),nsamp
cc      if(ill.eq.2) write(4,2016) (pbc(it),it=ip1,ip2),pbct
cc      write(4,2001) (mxx(5),i=1,ncpr)
cc      write(4,2002)
      if(ip2.eq.kx1) go to 440
      ip1=ip2+1
      ip2=kx1
      go to 445
  440 continue
cc      call aicp1(lk5,itemx,lc,aic,iw(i13),iw(i16),iw(i17),data,icon2,
cc     2           w(j12),w(j13),ikr,ikn,ikf,n,nsamp,iw(i20),nnn,in)
      call aicp10(lk5,itemx,lc,aic,data,icon2,ikr,ikn,n,nsamp,ier)
      if( ier(1).ne.0 ) return

cc      write(4,2045) aic
      aic1(ncount)=aic
cc      write(4,2034)
      do 550 i=1,lk6
      if(i.ne.1) ij=lc(i-1)
      if(i.eq.1) ij=1
      ijj=imm(ij)
      iii=idt(i)
      if(itype(ij).ne.1) go to 555
      itm=item(ijj)
      nn2=1
      nn1=0
      do 556 ii=2,itm
      if(icon(ij,ii-1).eq.icon2(ij,ii)) nn2=nn2+1
      if(icon(ij,ii-1).eq.icon2(ij,ii)) go to 556
      nn1=nn1+1
      ias2(ij,nn1)=nn2
      nn2=1
  556 continue
      nn1=nn1+1
      ias2(ij,nn1)=nn2
  555 continue
      is=1
         nint=0
      if(itx.ne.1.or.xx(ijj).eq.0.) is=ias2(ij,1)
      do 560 ii=1,iii
         nint=nint+1
      ie=is+ias2(ij,ii)
      if(itx.ne.1.or.xx(ijj).eq.0.) ie=is+ias2(ij,ii+1)
      if(itx.ne.1.and.item(ijj).le.is) is=item(ijj)
      if(itx.ne.1.and.ii.eq.iii) is=item(ijj)
      if(item(ijj)+1.le.ie) ie=item(ijj)+1
cc      if(ii.ne.1) go to 567
cc      if(i.ne.1) go to 565
cc      if(itx.ne.1.or.xx(ijj).eq.0.) write(4,2020) i,(c(j,1),j=1,20),
cc     &                                              ii,iaa(ijj,is)
cc      if(itx.eq.1.and.xx(ijj).ne.0.) write(4,2021) i,(c(j,i),j=1,20),
cc     &                                  ii,ab(ijj,is),ab(ijj,ie)
cc      go to 566
cc  565 if(itx.ne.1.or.xx(ijj).eq.0.) write(4,2020) i,
cc     &                            (cb(j,i-1),j=1,20),ii,iaa(ijj,is)
cc      if(itx.eq.1.and.xx(ijj).ne.0.) write(4,2021) i,
cc     &                (cb(j,i-1),j=1,20),ii,ab(ijj,is),ab(ijj,ie)
cc      go to 566
cc  567 continue
cc      if(itx.ne.1.or.xx(ijj).eq.0.) write(4,2040) ii,iaa(ijj,is)
cc      if(itx.eq.1.and.xx(ijj).ne.0.) write(4,2041) ii,
cc     &                                     ab(ijj,is),ab(ijj,ie)
cxx  566 is=ie
  566 continue
c <<<
      if(itx.ne.1.or.xx(ijj).eq.0.) iabse(i,nint,ncount)=is
      if(itx.eq.1.and. xx(ijj).ne.0.) then
         if (nint .eq. 1) then
            iabse(i,nint,ncount)=is
            nint=nint+1
         end if
         iabse(i,nint,ncount)=ie
      end if
      is=ie
c >>>
  560 continue
cc      write(4,2002)
  550 continue
cc      if(lk6.ne.1.and.izu.eq.1) call pr1(ib,idt,iby,kk,ikk,n,lk6)
cc      if(lk6.eq.1.or.izu.ne.1) go to 751
cc      write(4,2034)
cc      do 750 i=1,lk6
cc      if(i.ne.1) ij=lc(i-1)
cc      if(i.eq.1) ij=1
cc      ijj=imm(ij)
cc      iii=idt(i)
cc      is=1
cc      if(itx.ne.1.or.xx(ijj).eq.0.) is=ias2(ij,1)
cc      do 760 ii=1,iii
cc      ie=is+ias2(ij,ii+1)
cc      if(itx.ne.1.or.xx(ijj).eq.0.) ie=is+ias2(ij,ii+1)
cc      if(itx.ne.1.and.item(ijj).le.is) is=item(ijj)
cc      if(itx.ne.1.and.ii.eq.iii) is=item(ijj)
cc      if(item(ijj)+1.le.ie) ie=item(ijj)+1
cc      if(ii.ne.1) go to 767
cc      if(i.ne.1) go to 765
cc      if(itx.ne.1.or.xx(ijj).eq.0.) write(4,2020) i,(c(j,1),j=1,20),
cc     &                                              ii,iaa(ijj,is)
cc      if(itx.eq.1.and.xx(ijj).ne.0.) write(4,2021) i,(c(j,i),j=1,20),
cc     &                                  ii,ab(ijj,is),ab(ijj,ie)
cc      go to 766
cc  765 if(itx.ne.1.or.xx(ijj).eq.0.) write(4,2020) i,
cc     &                            (cb(j,i-1),j=1,20),ii,iaa(ijj,is)
cc      if(itx.eq.1.and.xx(ijj).ne.0.) write(4,2021) i,
cc     &                (cb(j,i-1),j=1,20),ii,ab(ijj,is),ab(ijj,ie)
cc      go to 766
cc  767 continue
cc      if(itx.ne.1.or.xx(ijj).eq.0.) write(4,2040) ii,iaa(ijj,is)
cc      if(itx.eq.1.and.xx(ijj).ne.0.) write(4,2041) ii,
cc     &                                     ab(ijj,is),ab(ijj,ie)
cc  766 continue
cc      is=ie
cc  760 continue
cc      write(4,2002)
cc  750 continue
cc  751 continue
cc      write(4,2043)
      if(icl.eq.0.or.ncount.gt.icl) go to 580
cc      if(ncount.lt.10) write(4,2044) ncount
cc      if(ncount.ge.10) write(4,2046) ncount
      do 620 i=ncount,icl
      lc1=imm(1)
      if(icls(2,i).eq.lc1) go to 630
  620 continue
      go to 580
  630 ncount=ncount+1
      ii=i
      lk5=icls(1,ii)-1
cx          if((lk5+2).gt.10) go to 588
          if((lk5+2).gt.max(10,n+1)) go to 588
      do 640 i=1,lk5
      lk9=icls(i+2,ii)
      do 650 j=1,n
      if(lk9.eq.imm(j)) go to 660
  650 continue
cc      write(4,*) 'err'
      ier(1)=650
      go to 580
  660 lca(i)=j
  640 continue
      go to 600
  580 continue
cc      write(4,2017)
      return
cc  570 write(4,2022) ialim,imax1,jalim,jmax1
  570 continue
         ier(1)=2022
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
      return
c
  573 continue
         ier(1)=2003
         ier(2)=icaa(i)
         return
  574 continue
         ier(1)=2003
         ier(2)=lk4
         return
c
cc  575 write(4,2048) ikkk,ikk
  575 continue
         ier(1)=2048
         ier(2)=ikk
         return
  585 continue
cddd         ier(1)=2048
         ier(2)=kkj
         return
  588 continue
         ier(1)=2588
         ier(2)=lk5+2
c
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
      return
 2001 format(15x,100a1)
 2002 format(' ')
 2003 format(' ')
 2004 format(15x,'response variable   (',20a1,')')
 2005 format('+',53('-'),'+'/
     1       'i contingency table constructed by the best subset of i'/
     2       'i explanatory variables                               i'/
     3       '+',53('-'),'+')
 2006 format(2a1,i1,2a1,20a1)
 2007 format(1h+,28x,3(2a1,i1,2a1,20a1,2x)/
     1       (28x,3(2a1,i1,2a1,20a1,2x)))
 2008 format(1x,10a1,5x,'response variable')
 2009 format(1h+,27x,' x(1) ')
 2010 format(15x,20i5)
 2011 format(24x,'    codes        class bound',
     1'aries')
 2012 format(1x,10i1)
 2013 format(1h+,30x,20i5)
 2014 format(1h+,30x,16f6.1)
 2015 format(' total',10x,20i5)
 2016 format(' total',10x,16f6.1)
 2017 format(' ')
 2018 format(i3,1x,(20a1,1x),i4,2x,f10.2,2x,f10.2,2x,f10.5)
 2019 format((4x,(20a1,1x)))
 2020 format(' x(',i1,'):',20a1,1x,i5,' : ',i5)
 2021 format(' x(',i1,'):',20a1,1x,i5,' : ',d12.5,' - ',d12.5)
 2022 format(' ia or a dimension over ',4i10)
 2023 format(72('-'))
 2024 format('number of explanatory variables = ',i5)
 2026 format(67('-')/'           ', 11x,'number of   ',
     1 2x,'     '/
     2 4x,'explanatory',7x,'categories    ',10x,'difference'/
     3 4x,'variables',9x,'of exp. var.',2x,'a i c',5x,'of aic',7x,
     4 'weight'/67('-'))
 2027 format(' response variable  : ',20a1/)
 2028 format(72('=')//'+',62('-'),'+'/
     1        'i aic''s of the models with k explanatory varia',
     2        'bles (k=1,2,...) i'/'+',62('-'),'+')
 2029 format('+',45('-'),'+'/
     1       'i summary of subsets of explanatory variables i'/
     2        '+', 45('-'),'+')
 2030 format(///10x,'---------------        ----------------'
     1         /10x,'---------------  note  ----------------'
     2         /10x,'the simple minded choice of the explanatory ',
     3                  'variables with'
     3         /10x,'the minimum of aic may not ',
     4                  'produce the best result. '
     5         /10x,'rather it is recommended to observe the ',
     6                  'general behavior of aic''s.'
     7         /10x,'---------------  note  ----------------'
     8         /10x,'---------------        ----------------')
 2031 format(1h ,15x,20a1)
 2032 format(1h+,10x,' - - -')
 2033 format(15x,16i6)
 2034 format(' < note > ')
 2035 format(/'<5',i1,'>')
 2036 format(/'<5>')
 2037 format(72('=')//'<6>')
 2040 format(27x,i5,' : ',i5)
 2041 format(27x,i5,' : ',d12.5,' - ',d12.5)
 2042 format(67('-'))
 2043 format(72('+'))
 2044 format(/'<7',i1,'>')
 2045 format(/'a i c =',f10.2/)
 2046 format(/'<7',i2,'>')
 2047 format('+---------------------------------------+'/
     1       'i the output of the additional analysis i'/
     2       '+---------------------------------------+')
 2048 format(' ikkk over  ',2i10)
      end
cc      subroutine aicp(k,idf,lc,lk2,aicc,ld,le,lb,ni,l1,l2,ly,lp,c,
cc     1                cb,lz,ca,ica,nni,data,icon,am,acmm,bic,ac,aa,
cc     2                lk7,ikr,ikn,ikf,n,nsamp,iw,nnn,in)
      subroutine aicp0(k,idf,lc,lk2,aicc,ca,ica,nni,data,icon,aa,lk7,
     1                 ikr,jkr,ikn,n,nsamp,eps01,ier)
c
c     this subroutine provides various combinations of explanatory
c     variables to compute the corresponding aic's.   this is almost
c     the same as the subroutine 'aicm'.
c
      implicit real*8(a-h,o-z)
cc      integer *2 lc,ld,le,lb,idf,l1,l2,lp,ly,c,cb,lz,ni,ca,ica,iw,
cc     1           nni,data,icon,nnmm
      integer c,cb,ca,data
      dimension lc(n),ld(n,n),am(n),acmm(n),le(n,n),bic(n),
     1          lb(n),ni(n),idf(n),l1(n),l2(n),ly(2,n),lp(2),ac(n),
     2          lz(n),aa(ikr),icon(n,ikn),
     3          ca(ikr,jkr),ica(ikr),nni(ikr),data(nsamp,n),ier(2)
cc      common ialim,imax,jalim,jmax,imaxx,jmaxx
c
      lk7=0
      do 10 i=1,k
      l1(i)=i
   10 continue
      do 20 l=1,n
      lz(l)=0
   20 continue
cc      i11=1+n
cc      i12=i11+n
cc      i13=i12+2
cc      i14=i13+ikf*2
cc      i15=i14+ikf
cc      i16=i15+ikf
cc      imax1=imax+i16
cc      jmax1=jmax
cc      if(imax1.gt.ialim) go to 330
cc      if(imax1.gt.imaxx) imaxx=imax1
cc      imax=imax1
      lk=0
      kk=k
      lk4=1
   30 continue
      do 31 i=1,2
      do 31 j=1,n
   31 ly(i,j)=0
      ly(1,1)=1
      lp(1)=1
      if(lk.ne.0) go to 45
      lk3=1
c <<<
      ikf=idf(1)
c >>>
cc      call aicsub (idf,lk3,ac(1),ni(1),bic(1),iw(1),iw(i11),iw(i12),
cc     1             iw(i13),ly,lp,data,icon,nsamp,n,ikf,ikn,
cc     2             iw(i14),iw(i15),in)
      call aicsub0 (idf,lk3,ac(1),ni(1),bic(1),ly,lp,data,icon,nsamp,n,
     1             ikf,ikn,ier)
      if( ier(1).ne.0 ) return
c-----    modified by M.I.
      if(ac(1) .gt. aa(lk4)) go to 45
c-----
      aa(lk4)=ac(1)
      nni(lk4)=ni(1)
      ca(lk4,1)=l1(1)
      ica(lk4)=lk+1
   45 continue
      do 60 i=2,kk
      lk4=lk4+1
      if(lk4.gt.ikr) go to 340
      lk3=2
      ly(1,2)=l1(i)
      ly(2,1)=l1(i)
      if(lk.eq.0) go to 70
      do 80 l=1,lk
      ly(2,l+1)=lb(l)
      ly(1,l+2)=lb(l)
   80 continue
   70 continue
      lp(1)=lk+2
      lp(2)=lk+1
c <<<
      ikf=1
      lpp=lp(1)
      do 71 ijk=1,lpp
      lyy=ly(1,ijk)
      ikf=ikf*idf(lyy)
   71 continue
c >>>
cc      call aicsub (idf,lk3,ac(i),ni(i),bic(i),iw(1),iw(i11),iw(i12),
cc     1             iw(i13),ly,lp,data,icon,nsamp,n,ikf,ikn,
cc     2             iw(i14),iw(i15),in)
      call aicsub0 (idf,lk3,ac(i),ni(i),bic(i),ly,lp,data,icon,nsamp,n,
     1             ikf,ikn,ier)
      if( ier(1).ne.0 ) return
c-----    modified by M.I.
      if(ac(i) .gt. aa(lk4)) go to 65
c-----

      aa(lk4)=ac(i)
      nni(lk4)=ni(i)
      ca(lk4,1)=l1(1)
      if(lk.eq.0) go to 66
         if((lk+1).gt.jkr) go to 900
      do 67 l=1,lk
   67 ca(lk4,l+1)=lb(l)
         if((lk+2).gt.jkr) go to 900
   66 ca(lk4,lk+2)=l1(i)
      ica(lk4)=lk+2
   65 continue
   60 continue
      
      amin2=10.**10
      amin=10.**10
c <<<
      ikm=2
c >>>
      do 90 i=2,kk
      if(amin2.lt.bic(i)) go to 140
      amin2=bic(i)
  140 continue
ccxx      if(amin.lt.ac(i)) go to 90
         ddd=dabs(amin-ac(i))
         if(ddd.ne.0) ddd=ddd/max(dabs(amin), dabs(ac(i)))
      if(amin.lt.ac(i) .and. ddd.gt.eps01) go to 90

      amin=ac(i)
      ikm=i
   90 continue
      jj=ikm
      lk=lk+1
      lb(lk)=l1(jj)
      ii=0
      acmm(lk)=amin
      do 150 i=1,lk
      le(lk,i)=lb(i)
  150 continue
      do 160 i=1,kk
      if(jj.eq.i) go to 160
      ii=ii+1
      l2(ii)=l1(i)
  160 continue
      kk=ii
      if(kk.eq.1) go to 170
      do 180 i=1,kk
      l1(i)=l2(i)
  180 continue
      if(lk.lt.3) go to 30
      lkk=lk-2
c <<<
      lk1=lk-1
c >>>
      do 190 i=1,lkk
      jj=0
      do 200 j=1,lk
      if(i.eq.j) go to 200
      jj=jj+1
      lc(jj)=lb(j)
      ld(i,jj)=lb(j)
  200 continue
cc      lk1=lk-1
      do 210 l=1,lk1
      ly(1,l+1)=lc(l)
      ly(2,l)=lc(l)
  210 continue
      lp(1)=lk
      lp(2)=lk1
      lk3=2
c <<<
      ikf=1
      lpp=lp(1)
      do 211 ijk=1,lpp
      lyy=ly(1,ijk)
      ikf=ikf*idf(lyy)
  211 continue
c >>>
cc      call aicsub (idf,lk3,am(i),nnmm,ax,iw(1),iw(i11),iw(i12),
cc     1             iw(i13),ly,lp,data,icon,nsamp,n,ikf,ikn,
cc     2             iw(i14),iw(i15),in)
      call aicsub0 (idf,lk3,am(i),nnmm,ax,ly,lp,data,icon,nsamp,n,ikf,
     1             ikn,ier)
      if( ier(1).ne.0 ) return

      lk4=lk4+1
      if(lk4.gt.ikr) go to 340
c-----   modified by M.I.
      if(am(i) .gt. aa(lk4)) go to 190
c-----
      aa(lk4)=am(i)
      nni(lk4)=nnmm
      ica(lk4)=lk1+1
         if((lk1+1).gt.jkr) go to 910
      do 215 j=1,lk1
  215 ca(lk4,j+1)=lc(j)
      ca(lk4,1)=l1(1)
  190 continue
      aminx=amin
      do 220 i=1,lkk
      if(aminx.le.am(i)) go to 220
      aminx=am(i)
      jj=i
  220 continue
      if(aminx.ne.amin.and.aminx.le.acmm(lk1)) go to 230
      if(acmm(lk1).le.acmm(lk)) go to 170
      go to 30
  230 do 240 i=1,lk1
      lb(i)=ld(jj,i)
      le(lk1,i)=ld(jj,i)
  240 continue
      ii=0
      do 250 i=1,k
      do 260 j=1,lk1
      if(lb(j).eq.i) go to 250
  260 continue
      ii=ii+1
      l1(ii)=i
  250 continue
      kk=ii
      acmm(lk1)=aminx
      lk=lk-1
      go to 30
  170 continue
      if(lk7.le.lk4) lk7=lk4
      jj=ikm
      do 270 ij=1,lk
      aminx=10.**10
      do 280 i=1,lk
      if(lz(i).ne.0) go to 280
      if(aminx.le.acmm(i)) go to 280
         daic=dabs(aminx-acmm(i))
         if(daic.ne.0) daic=daic/max(dabs(aminx), dabs(acmm(i)))
         if(daic.lt.eps01) go to 280

      aminx=acmm(i)
      lk3=i
  280 continue
      lz(lk3)=lk3
      if(ij.ne.1) go to 290
      lk2=lk3
cxx      do 300 ix=1,10
cxx      cb(ix,1)= c(ix,1)
cxx  300 continue
      do 310 i=1,lk2
      lc(i)=le(lk2,i)
  310 continue
  290 continue
cxx      do 320 i=1,lk3
cxx      lx=le(lk3,i)
cxx      iy=i+1
cxx      do 320 ix=1,10
cxx      cb(ix,iy )=c(ix,lx)
cxx  320 continue
      if(ij.eq.1)aicc=acmm(lk3)
  270 continue
      return
cc  330 write(4,2001) ialim,imax1,jalim,jmax1
  330 ier(1)=2001
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
      return
cc  340 write(4,2002) lk4,ikr
  340 continue
         ier(1)=2002
         ier(2)=lk4
         return
  900 continue
         ier(1)=2003
         ier(2)=lk
         return
  910 continue
         ier(1)=2003
         ier(2)=lk1
c
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
      return
 2001 format(' ia or a dimension over ',4i10)
 2002 format(' ikr over ',2i10)
      end
      subroutine bun0(icon,ias1,jj,aic,itemt,ii,itemx,lk5,ca,ica,
cc     1               nni,c,aa,n,iw,w,data,nsamp,ikr,ikn,ikf,lk77,nnn,in)
     1               nni,aa,n,data,nsamp,ikr,jkr,ikn,lk77,eps01,ier)
c
c     this subroutine re-divides the pooled categories of an
c     explanatory variable in a multidimensional table to search for
c     the maice.
c
      implicit real*8(a-h,o-z)
cc      integer*2 itemx,ias1,icon,lc,itemt,c,ca,ica,nni,data,iw
      integer ca,data
      dimension icon(n,20),ias1(2,n,20),itemx(n),data(nsamp,n),lc(n),
cc     1          c(10,n),ca(ikr,10),ica(ikr),nni(ikr),aa(ikr),itemt(2,n),
     1           ca(ikr,jkr),ica(ikr),nni(ikr),aa(ikr),itemt(2,n)
      dimension ier(2)
cc     3          iw(nnn),w(1)
cc      common ialim,imax,jalim,jmax,imaxx,jmaxx
cc      i11=1+n*n
cc      i12=i11+n*n
cc      i13=i12+n
cc      i14=i13+n
cc      i15=i14+n
cc      i16=i15+n
cc      i17=i16+2*n
cc      i18=i17+2
cc      i19=i18+10*n
cc      i20=i19+n
cc      j11=1+n
cc      j12=j11+n
cc      j13=j12+n
cc      j14=j13+n
cc      imax1=imax+i20
cc      jmax1=jmax+j14
cc      if(imax1.gt.ialim.or.jmax1.gt.jalim) go to 120
c
      aic=1.d10
      ii=0
      do 10 i=2,n
      if(i.eq.jj) go to 10
      if(ias1(1,i,1).eq.1) go to 10
      do 30 k=1,n
      kx=1
      if(i.eq.k) kx=2
      itemx(k)=itemt(kx,k)
      k1=0
      itemm=itemx(k)
      do 40 kk=1,itemm
      iax=ias1(kx,k,kk+1)
      do 50 kkk=1,iax
      k1=k1+1
      icon(k,k1)=kk
   50 continue
   40 continue
   30 continue
cc      imax=imax1
cc      jmax=jmax1
cc      call aicp(n,itemx,lc,lk2,aic1,iw(1),iw(i11),iw(i12),iw(i13),
cc     1          iw(i14),iw(i15),iw(i16),iw(i17),c,iw(i18),iw(i19),
cc     2          ca,ica,nni,data,icon,w(1),w(j11),w(j12),w(j13),aa,
cc     3          lk7,ikr,ikn,ikf,n,nsamp,iw(i20),nnn,in)
      call aicp0(n,itemx,lc,lk2,aic1,ca,ica,nni,data,icon,aa,lk7,ikr,
     1          jkr,ikn,n,nsamp,eps01,ier)
      if( ier(1).ne.0 ) return

      if(aic.le.aic1) go to 10
         daic=dabs(aic-aic1)
         if(daic.ne.0) daic=daic/max(dabs(aic), dabs(aic1))
         if(daic.lt.eps01) go to 10

      aic=aic1
      ii=i
      lk5=lk2
      lk77=lk7
   10 continue
      if(ii.eq.0) go to 20
      do 130 k=1,n
      kx=1
      if(ii.eq.k) kx=2
      itemx(k)=itemt(kx,k)
      k1=0
      itemm=itemx(k)
      do 140 kk=1,itemm
      iax=ias1(kx,k,kk+1)
      do 150 kkk=1,iax
      k1=k1+1
      icon(k,k1)=kk
  150 continue
  140 continue
  130 continue
   20 continue
      return
cc  120 write(4,2001) ialim,imax1,jalim,jmax1
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
cc 2001 format(' ia or a dimension over ',4i10)
      end
cc      subroutine aicsub(idf,l,ac,ni,bc,idt,idd,id1,ia,ly,lp,data,icon,
cc     1                  nsamp,n,ikf,ikn,tc,tr,in)
      subroutine aicsub0(idf,l,ac,ni,bc,ly,lp,data,icon,nsamp,n,ikf,ikn,
     1 ier)
c
c     this subroutine computes aic values for a (multidimensional)
c     contingency table.
c
      implicit real*8 (a-h,o-z)
cc      integer*2 ia,idf,lp,ly,ni,data,idt,idd,id1,icon,tc,tr
      integer data,tc,tr
      dimension idt(n),idd(n),id1(2),ia(ikf,2),idf(n),lp(2),ly(2,n),
     1          data(nsamp,n),icon(n,ikn),tc(ikf),tr(ikf),ier(2)
      t=nsamp
      expo=exp(-1.)
      do 10 i=1,ikf
      tc(i)=0
      tr(i)=0
      do 10 j=1,2
   10 ia(i,j)=0
      iaa=1
      lpp=lp(1)
      do 15 i=1,lpp
      lyy=ly(1,i)
      iaa=iaa*idf(lyy)
   15 continue
      if(iaa.gt.ikf) go to 120
      do 20 i=1,nsamp
      do 30 ii=1,l
      lpp=lp(ii)
      if(lpp.eq.0) go to 30
      do 40 iii=1,lpp
      lyy=ly(ii,iii)
      idf1=idf(lyy)
      idt(iii)=idf1
      id9=data(i,lyy)
      idd(iii)=id9
   40 continue
      do 41 iii=1,lpp
      lyy=ly(ii,iii)
      idf1=idf(lyy)
      idt(iii)=idf1
      id9=data(i,lyy)
      id=icon(lyy,id9)
      idd(iii)=id
   41 continue
      iaa=idd(1)
      if(lpp.eq.1) go to 50
      ii2=1
      do 60 iii=2,lpp
      ii2=ii2*idt(iii-1)
      iaa=iaa+ii2*(idd(iii)-1)
   60 continue
   50 continue
      ia(iaa,ii)=ia(iaa,ii)+1
      if(ii.eq.2) tr(iaa)=tr(iaa)+1
   30 continue
      ii3=ly(1,1)
      id=data(i,ii3)
      tc(id)=tc(id)+1
   20 continue
      do 70 i=1,l
      id1(i)=1
      lpp=lp(i)
      if(lpp.eq.0) go to 80
      do 70 j=1,lpp
      lyy=ly(i,j)
      id1(i)=id1(i)*idf(lyy)
   70 continue
   80 continue
      ic0=0
      idf1=idf(1)
      do 21 i=1,idf1
      if(tc(i).eq.0) ic0=ic0+1
   21 continue
      ir0=0
      if(l.eq.1) go to 23
      idf1=id1(2)
      do 22 i=1,idf1
      if(tr(i).eq.0) ir0=ir0+1
   22 continue
   23 continue
      ip=0
      c=0.
c <<<
      t1=0.
      i61=0
c >>>
      do 100 j=1,l
      kk=id1(j)
      if(j.eq.2) go to 95
      t1=0.
      kkk=idf(1)
      kk1=kk/kkk
      do 90 m=1,kk1
      do 91 k=1,kkk
      ii=(m-1)*kkk+k
      if(tc(k).eq.0) go to 91
      if(l.ne.1.and.tr(m).eq.0) go to 91
      if(ia(ii,1).eq.0) t1=t1+expo
   91 continue
   90 continue
   95 continue
      a=0.
      if(kk.eq.1) go to 100
      do 110 k=1,kk
      aaa=ia(k,j)
      if(aaa.eq.0.) aaa=expo
      if(j.eq.1) go to 114
      if(tr(k).eq.0) go to 110
      kkk=idf(1)
      kk1=0
      do 125 m=1,kkk
      if(tc(m).eq.0) go to 125
      kk2=(k-1)*kkk+m
      if(ia(kk2,1).eq.0) kk1=kk1+1
  125 continue
      aaa=kk1*expo+ia(k,j)
      go to 115
  114 continue
      kk2=(k-1)/idf(1)+1
      if(l.ne.1.and.tr(kk2).eq.0) go to 110
      idf1=idf(1)
      kk2=mod(k,idf1)
      if(kk2.eq.0) kk2=idf(1)
      if(tc(kk2).eq.0) go to 110
  115 continue
      a=aaa*log(aaa/(t+t1))+a
  110 continue
      if(j.eq.1.and.l.eq.1) i61=(idf(1)-ic0)-1
      if(j.eq.1.and.l.ne.1) i61=(idf(1)-ic0)*(id1(2)-ir0)-1
      if(j.eq.2) i61=(kk-ir0)-1
      if(j.eq.2) a=-a
      if(j.eq.2) i61=-i61
      c=c+a
      ip= i61+ip
  100 continue
      a=0.
      kkk=idf(1)
      kk=id1(1)/kkk
      do 130 i=1,kkk
      if(tc(i).eq.0) go to 130
      aaa=0.
      do 140 j=1,kk
      if(l.ne.1.and.tr(j).eq.0) go to 140
      ii=(j-1)*kkk+i
      if(ia(ii,1).eq.0) aaa=aaa+expo
      aaa=aaa+ia(ii,1)
  140 continue
      a=aaa*log(aaa/(t+t1))+a
  130 continue
      c=c-a
      i61=kkk-ic0-1
      ip=ip-i61
      ni=id1(2)
      if(l.eq.1) ni=0
      ac  = -2*(c-ip)
      bc = ac
      return
cc  120 write(4,2001) iaa,ikf
  120 continue
         ier(1)=2001
         ier(2)=iaa
cc      close(1)
cc      close(4)
cc      close(in)
 2001 format(' dimension over ia ' ,2i10)
cc      stop 10
      return
      end
cc      subroutine aicp1(k,idf,lc,aicc,ni,ly,lp,data,icon,bic,ac,
cc     1                ikr,ikn,ikf,n,nsamp,iw,nnn,in)
      subroutine aicp10(k,idf,lc,aicc,data,icon,ikr,ikn,n,nsamp,ier)
c
c     this subroutine computes the value of aic for the output <7>.
c
      implicit real*8(a-h,o-z)
cc      integer *2 lc,idf,lp,ly,ni,iw,data,icon
      integer data
      dimension lc(n),bic(n),ni(n),idf(n),ly(2,n),lp(2),ac(n),
cc     1          icon(n,ikn),iw(nnn),data(nsamp,n)
     1          icon(n,ikn),data(nsamp,n),ier(2)
cc      common ialim,imax,jalim,jmax,imaxx,jmaxx
cc      i11=1+n
cc      i12=i11+n
cc      i13=i12+2
cc      i14=i13+ikf*2
cc      i15=i14+ikf
cc      i16=i15+ikf
cc      imax1=imax+i16
cc      jmax1=jmax
cc      if(imax1.gt.ialim) go to 330
cc      if(imax1.gt.imaxx) imaxx=imax1
cc      imax=imax1
cc      lk4=1
cc   30 continue
      do 31 i=1,2
      do 31 j=1,n
   31 ly(i,j)=0
      ly(1,1)=1
      lp(1)=1
      lk3=2
      do 80 l=1,k
      ly(2,l)=lc(l)
      ly(1,l+1)=lc(l)
   80 continue
      lp(1)=k+1
      lp(2)=k
      i=1
cc      call aicsub (idf,lk3,ac(i),ni(i),bic(i),iw(1),iw(i11),iw(i12),
cc     1             iw(i13),ly,lp,data,icon,nsamp,n,ikf,ikn,
cc     2             iw(i14),iw(i15),in)
c <<<
      ikf=1
      lpp=lp(1)
      do 81 ii=1,lpp
      lyy=ly(1,ii)
      ikf=ikf*idf(lyy)
   81 continue
c >>>
      call aicsub0 (idf,lk3,ac(i),ni(i),bic(i),ly,lp,data,icon,nsamp,n,
     1             ikf,ikn,ier)
      if( ier(1).ne.0 ) return

      aicc=ac(i)
      return
cc  330 write(4,2001) ialim,imax1,jalim,jmax1
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
cc  340 write(4,2002) lk4,ikr
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
cc 2001 format(' ia or a dimension over ',4i10)
cc 2002 format(' ikr over ',2i10)
      end
c
c
c--------------   common subroutine   --------------
c
cc      subroutine ac1p (a,ii,total,jx,nd,lc,lcc,lcy,lj,knn,lb,a1,totl,
cc     1                 tttr,aa,am,kk5,kk2,iw,w,n11,n33,totalc)
      subroutine ac1p (a,ii,total,jx,lc,lcc,lcy,lj,knn,a1,totl,tttr,aa,
     1                  am,kk5,kk2,n11,n33,totalc,eps01)
c
c     this subroutine searches for maice within the possible ways of
c     categorization of a single explanatory variable in cace
c     'ity(k)=1'.
c
      implicit real*8(a-h,o-z)
cc      integer*2 nd,lc,lcc,knn,lb,total,lcy,lj,totalc,iw
      integer total,totalc
      dimension a(n11,n33),a1(n11,n33),total(n33),totl(n33),tttr(n33),
     1          nd(n33),lc(n33),lcc(n33),lcy(10,n33),aa(n33),lj(10),
cc     2          knn(n33),lb(n33),iw(1),w(1),totalc(n11)
     2          knn(n33),lb(n33),totalc(n11)
cc      common ialim,imax,jalim,jmax,imaxx,jmaxx
c
c     substituting exp(-1.0) for zero frequency
c
      expo=exp(-1.)
      al=0.
      jj=jx
      samp=0.0
      tsmp=0.0
      do 15 j=1,jj
      tttr(j)=0.
      lc(j)=1
      lcy(1,j)=1
      if(total(j).eq.0) go to 15
      do 10 i=1,ii
      if(totalc(i).eq.0) go to 10
      if(a(i,j).le.0.0) tsmp=tsmp+expo
      if(a(i,j).le.0.) tttr(j)=tttr(j)+expo
      samp=samp+a(i,j)
   10 continue
   15 continue
c
c     computation of aic for the original two-way table
c
      do 21 i=1,ii
      if(totalc(i).eq.0) go to 21
      tt=0.0
      do 25 j=1,jj
      if(total(j).eq.0) go to 25
      if(a(i,j).le.0.0) tt=tt+expo
      tt=tt+a(i,j)
   25 continue
      do 20 j=1,jj
      if(total(j).eq.0) go to 20
      aaa=a(i,j)
      if(aaa.le.0.) aaa=expo
      al=al+aaa*log(aaa/(tt*(total(j)+tttr(j)))*(samp+tsmp))
   20 continue
   21 continue
      ii0=ii
      do 22 i=1,ii
   22 if(totalc(i).eq.0) ii0=ii0-1
      jj0=jj
      do 23 j=1,jj
   23 if(total(j).eq.0) jj0=jj0-1
      lj(1)=1
c     in=ii*jj-jj-ii+1
      in=ii0*jj0-jj0-ii0+1
      aa(1)=-2.*(al-in)
      nd(1)=in
      kx=jj
      knn(1)=kx
      km=2
c <<<
      k2=2
c >>>
   30 continue
c
c     pooling of adjacent categories of an explanatory variable,
c     computation of the corresponding aic and the search for maice:
c       (the number of categories pooled increases in stepwise .)
c
      if(kx.le.1) go to 40
      aicmi=10.**10
      if(kx.eq.2) go to 50
      do 60 i=2,kx
      do 70 j=1,kx
      if(j.lt.i) lcc(j)=lc(j)
      if(j.eq.i) lcc(j-1)=lcc(j-1)+lc(j)
      if(j.gt.i) lcc(j-1)=lc(j)
   70 continue
      jj=kx-1
      k1=jj
      j1=0
      do 80 j2=1,jj
      totl(j2)=0.
   80 continue
      do 90 j2=1,jj
      j3=lcc(j2)
      do 100 m=1,ii
      a1(m,j2)=0.
  100 continue
      do 90 j4=1,j3
      j1=j1+1
      totl(j2)=totl(j2)+total(j1)
      do 90 m=1,ii
      a1(m,j2)=a1(m,j2)+a(m,j1)
   90 continue
      tsmp=0.0
      do 111 j2=1,jj
      if(totl(j2).eq.0.) go to 111
      tttr(j2)=0.
      do 110 m=1,ii
      if(totalc(m).eq.0) go to 110
      if(a1(m,j2).le.0.0) tsmp=tsmp+expo
      if(a1(m,j2).le.0.) tttr(j2)=tttr(j2)+expo
  110 continue
  111 continue
      al=0.
      do 120 m=1,ii
      if(totalc(m).eq.0) go to 120
      tt=0.0
      do 125 j2=1,jj
      if(totl(j2).eq.0.) go to 125
      if(a1(m,j2).le.0.0) tt=tt+expo
      tt=tt+a1(m,j2)
  125 continue
      do 130 j2=1,jj
      if(totl(j2).eq.0.) go to 130
      aaa=a1(m,j2)
      if(aaa.le.0.0) aaa=expo
      al=al+aaa*log(aaa/(tt*(totl(j2)+tttr(j2)))*(samp+tsmp))
  130 continue
  120 continue
      jj0=jj
      do 123 j2=1,jj
  123 if(totl(j2).eq.0) jj0=jj0-1
      in=ii0*jj0-jj0-ii0+1
      aic=-2.*(al-in)
ccxx      if(aicmi.lt.aic) go to 60
         daic=dabs(aicmi-aic)
         if(daic.ne.0) daic=daic/max(dabs(aicmi), dabs(aic))
      if(aicmi.lt.aic .and. daic.gt.eps01) go to 60

      iix=i-1
      aicmi=aic
      do 140 j=1,k1
      lb(j)=lcc(j)
  140 continue
   60 continue
      aicmm=aa(km-1)
      k2=k1
      if(lc (iix).le.3) go to 150
      do 160 j=1,k1
      lcc(j)=lb(j)
  160 continue
cc      i11=1+n33
cc      i12=i11+n33
cc      j11=1+n11*n33
cc      j12=j11+n33
cc      j13=j12+n33
cc      imax1=imaxx+i12
cc      jmax1=jmaxx+j13
cc      imax0=imaxx
cc      jmax0=jmaxx
cc      if(imax1.gt.ialim.or.jmax1.gt.jalim) go to 280
c
c     re-dividing the pooled categories and computing the
c     corresponding aic
c
cc      call yy(a,lcc,iix,aicmm,jj,total,ii,k1,iw(1),iw(i11),w(1),w(j11),
cc     1        w(j12),n11,n33,totalc)
      call yy(a,lcc,iix,aicmm,jj,total,ii,k1,n11,n33,totalc,eps01)
cc      imaxx=imax0
cc      jmaxx=jmax0
  150 continue
      if(km.eq.1) go to 50
      if(aa(km-1).le.aicmm) go to 50
         daic=dabs(aa(km-1)-aicmm)
         if(daic.ne.0) daic=daic/max(dabs(aa(km-1)), dabs(aicmm))
         if (daic.lt.eps01) go to 50

      kx=k1
      km=km-1
      k2=k1
      do 170 j=1,k1
      lb(j)=lcc(j)
  170 continue
   50 continue
      k1=k2
      kx=k1
      aa(km)=aicmi
      nd(km)=in
      knn(km)=k1
      if(km.le.10) go to 180
      do 190 i=2,10
      i1=i-1
      ljj=lj(i)
      k1=knn(ljj)
      do 200 j=1,k1
      lcy(i1,j)=lcy(i,j)
  200 continue
      lj(i1)=lj(i)
  190 continue
      k1=knn(km)
      do 210 j=1,k1
      lcy(10,j)=lb(j)
      lc(j)=lb(j)
  210 continue
      lj(10)=km
      go to 220
  180 continue
      do 230 i=1,k1
      lc(i)=lb(i)
      lcy(km,i)=lb(i)
  230 continue
      lj(km)=km
  220 continue
      if(km.le.2) go to 240
      if(aa(km).gt.aa(km-2).and.aa(km-1).gt.aa(km-2)) go to 250
  240 continue
      if(km.eq.n33) go to 250
      km=km+1
      go to 30
c
c     picking up the categorization with maice
c
   40 km=km-1
  250 continue
      am=10.**10
c <<<
      kk1=1
c >>>
      do 260 i=1,km
      if(am.lt.aa(i))go to 260
      am=aa(i)
      kk1=i
  260 continue
      kk5=knn(kk1)
      do 270 j=1,10
      if(kk1.eq.lj(j)) kk2=j
  270 continue
      return
cc  280 write(2,2001) ialim,imax,jalim,jmax
cc      close(1)
cc      close(4)
cc      close(in)
cc      stop 10
 2001 format(' ia or a dimension over ',4i10)
      end
cc      subroutine yy(a,lc,iix,aicmi,jj,total,ii,k1,lcc,lb,a1,totl,tttr,
cc     1              n11,n33,totalc)
      subroutine yy(a,lc,iix,aicmi,jj,total,ii,k1,n11,n33,totalc,eps01)
c
c     this subroutine divides the pooled categories of an explanatory
c     variable and computes the corresponding aic.
c
      implicit real*8(a-h,o-z)
cc      integer*2 lc,lcc,lb,total,totalc
      integer total,totalc
      dimension a(n11,n33),a1(n11,n33),total(n33),totl(n33),lc(n33),
     1          lcc(n33),lb(n33),tttr(n33),totalc(n11)
      amin=10.**10
c <<<
      k2=k1
c >>>
      ix=lc (iix)
      ix1=iix-1
      ix2=iix+1
      if(ix1.eq.0) go to 10
      do 20 i=1,ix1
      lcc(i)=lc(i)
   20 continue
   10 continue
      if(ix2.gt.k1) go to 35
      do 30 i=ix2,k1
      lcc(i+1)=lc(i)
   30 continue
   35 continue
      ix11=ix-1
      do 40 j=1,ix11
      lcc(iix)=j
      lcc(iix+1)=ix-j
      jj=k1+1
      j1=0
      do 50 j2=1,jj
      totl(j2)=0.
      j3=lcc(j2)
      do 60 m=1,ii
      a1(m,j2)=0.
   60 continue
      do 50 j4=1,j3
      j1=j1+1
      totl(j2)=totl(j2)+total(j1)
      do 50 m=1,ii
      a1(m,j2)=a1(m,j2)+a(m,j1)
   50 continue
      expo=exp(-1.)
      tsmp=0.0
      samp=0.0
      do 71 j2=1,jj
      if(totl(j2).eq.0.) go to 71
      tttr(j2)=0.0
      do 70 m=1,ii
      if(totalc(m).eq.0) go to 70
      if(a1(m,j2).le.0.0) tsmp=tsmp+expo
      if(a1(m,j2).le.0.) tttr(j2)=tttr(j2)+expo
      samp=samp+a1(m,j2)
   70 continue
   71 continue
      al=0.
      do 80 m=1,ii
      if(totalc(m).eq.0) go to 80
      tt=0.0
      do 85 j2=1,jj
      if(totl(j2).eq.0.) go to 85
      if(a1(m,j2).le.0.0) tt=tt+expo
      tt=tt+a1(m,j2)
   85 continue
      do 90 j2=1,jj
      if(totl(j2).eq.0.) go to 90
      aaa=a1(m,j2)
      if(aaa.le.0.0) aaa=expo
      al=al+aaa*log(aaa/(tt*(totl(j2)+tttr(j2)))*(samp+tsmp))
   90 continue
   80 continue
      ii0=ii
      do 81 m=1,ii
   81 if(totalc(m).eq.0) ii0=ii0-1
      jj0=jj
      do 82 j2=1,jj
   82 if(totl(j2).eq.0.) jj0=jj0-1
      in=ii0*jj0-jj0-ii0+1
      ac=-2.*(al-in)
      aic=ac
      if(amin.le.aic) go to 40
         daic=dabs(amin-aic)
         if(daic.ne.0) daic=daic/max(dabs(amin), dabs(aic))
         if(daic.lt.eps01) go to 40

      amin=aic
      k2=k1+1
      do 100 i=1,k2
      lb(i)=lcc(i)
  100 continue
   40 continue
      if(aicmi.le.amin) go to 110
         daic=dabs(aicmi-amin)
         if(daic.ne.0) daic=daic/max(dabs(aicmi), dabs(amin))
         if(daic.lt.eps01) go to 110

      k1=k2
      do 120 i=1,k2
      lc(i)=lb(i)
  120 continue
  110 continue
      return
      end
      subroutine eqck(caa,ikr,lk4,lk3,lk33,ijk)
cc      integer *2 caa,ca1,ca2,ca11,ca22
cc      dimension caa(ikr,10),ca1(10),ca2(10)
      integer caa,ca1,ca2,ca11,ca22
      dimension caa(ikr,lk4),ca1(lk4),ca2(lk4)
      ijk=0
      do 10 i=1,lk4
      ca1(i)=caa(lk3,i)
      ca2(i)=caa(lk33,i)
   10 continue
      do 20 i=1,lk4-1
      ca11=ca1(i)
      ca22=ca2(i)
      j1=i
      j2=i
      do 30 j=i+1,lk4
      if(ca11.lt.ca1(j)) go to 40
      ca11=ca1(j)
      j1=j
   40 if(ca22.lt.ca2(j)) go to 30
      ca22=ca2(j)
      j2=j
   30 continue
      ca1(j1)=ca1(i)
      ca1(i)=ca11
      ca2(j2)=ca2(i)
      ca2(i)=ca22
      if(ca11.ne.ca22) go to 50
   20 continue
      if(ca1(lk4).ne.ca2(lk4)) go to 50
      return
   50 ijk=1
      return
      end
