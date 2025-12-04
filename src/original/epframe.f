C-    *****************************************
C-
C-    ELASTIC PLASTIC ANALYSIS OF A PLANE FRAME
C-
C-                    PROGRAMMED BY HACKSOO LEE
C-                                FEB. 25, 1986
C-
C-    *****************************************
C-
C-	dimensionalize quantities
C-
      CHARACTER FILE7     
      CHARACTER*12 FNAME9,FNAME8,FNAME7
      DIMENSION A(600,900),ASAT(600,610),CORD(300,2),JTYPE(300,3) 
      DIMENSION SMA(300),AREA(300),PM(300),OLEN(300),VL(600),SF(600,2)
      DIMENSION SA(300),CSAT(900),SATX(900),ALG(600),CM(600),CT(300)
      DIMENSION MCON(300,2),CX(900)
      LOGICAL STORE
      STORE=.FALSE.
C-
C-	interactive user input - output
C-
    1 WRITE(*,1000)
 1000 FORMAT(/,'      OPTIONS TO RUN PROGRAM:',/,
     $           ' ','          0 - QUIT',/,
     $           ' ','          1 - READ DATA FROM SCREEN',/,
     $           ' ','          2 - READ DATA FROM DISK FILE',/)
      WRITE(*,'(A)')'      YOUR CHOICE ?'
      READ(*,'(I1)') IOPT
      IF (IOPT .EQ. 1) IDSK=5
      IF (IOPT .EQ. 2) IDSK=9
 1010 FORMAT(1X,/)
      IF (IOPT .EQ. 0) STOP '      EXECUTION ABORTED'
      IF (IOPT .EQ. 1) WRITE(*,'(A)')'      DO YOU LIKE TO STORE DATA F
     $ROM SCREEN IN A FILE (Y/N) ?'
      IF (IOPT .EQ. 1) READ(*,'(A)') FILE7
      IF ((IOPT.EQ.1) .AND. ((FILE7.EQ.'Y').OR.(FILE7.EQ.'y'))) 
     $ STORE=.TRUE.
      IF(STORE)WRITE(*,'(A)')'      ENTER FILE NAME FOR STORING DATA ?'
      IF (STORE) READ(*,'(A)') FNAME7
      IF (STORE) OPEN (7, FILE=FNAME7,STATUS='NEW')
      IF (IOPT .EQ. 2) WRITE(*,'(A)')'      DATA FILE NAME ?'
      IF (IOPT .EQ. 2) READ(*,'(A)') FNAME9
      IF (IOPT .GT. 2) GOTO 1
      IF (IOPT .EQ. 2) OPEN (IDSK, FILE=FNAME9)
      WRITE(*,'(A)') '      OUTPUT FILE NAME ?'
      READ(*,'(A)') FNAME8
      OPEN (8, FILE=FNAME8,STATUS='NEW')
      IF (IOPT .EQ. 1) WRITE (*,'(A)')'   *** START TO INPUT DATA ***'
      IF (IOPT .EQ. 1) WRITE (*,1010)
      IF (IOPT .EQ. 1) WRITE (*,'(A)')'      FRAME NUMBER ?'
      READ(IDSK,*) JFN
      IF (STORE) WRITE(7,'(I6)') JFN
      IF (IOPT .EQ. 1) WRITE(*,1020)
 1020 FORMAT(/,'      ENTER "JCT","NM", and "E" in one line',/,
     $      ' ','             1) JCT = NUMBER OF JUNCTIONS or JOINTS',/,
     $      ' ','             2) NM = NUMBER OF MEMBERS',/,
     $      ' ','             3) E = MODULUS OF ELASTICITY in ksi',/)
      READ(IDSK,*) JCT,NM,E
      IF (STORE) WRITE(7,'(2I6,F15.1)') JCT,NM,E
      WRITE(8,400) JFN
  400 FORMAT('%',5X,'ELASTIC PLASTIC ANALYSIS OF FRAME NO',I3/
     $       '%',5X,'---------------------------------------'/'%')
      IF (IOPT .EQ. 1) WRITE(*,1030)
 1030 FORMAT(/,'      ENTER "I","X","Y","DFX","DFY", and "DFZ" in one l
     $ine',/,'                              (repeated JCT times)',/,
     $      ' ','             1) I = JOINT NUMBER (from no.1)',/,
     $      ' ','             2) X = X-COORDINATE OF JOINT,inch',/,
     $      ' ','             3) Y = Y-COORDINATE OF JOINT,inch',/,
     $      ' ','             4) DFX = DEGREE OF FREEDOM IN THE X-DIRECT
     $ION (0 or 1)',/,
     $      ' ','             5) DFY = DEGREE OF FREEDOM IN THE Y-DIRECT
     $ION (0 or 1)',/,
     $      ' ','             6) DFZ = DEGREE OF FREEDOM IN ROTATION (0 
     $or 1)',/,
     $   ' ','                   ( 0 - RESTRAINED    1 - FREE )',/)
      DO 500 II=1,JCT
         READ(IDSK,*) I,(CORD(I,J),J=1,2),(JTYPE(I,J),J=1,3)
         IF (STORE) WRITE(7,1100) I,(CORD(I,J),J=1,2),(JTYPE(I,J),J=1,3)
  500 CONTINUE
 1100 FORMAT(' ',I6,2(F15.5,5X),3I5)
      IF (IOPT .EQ. 1) WRITE(*,1040)
 1040 FORMAT(//,'      ENTER "M","JT1","JT2","MI","A", and "MP" in one l
     $ine',/,'                               (repeated NM times)',/,
     $      ' ','             1) M = MEMBER NUMBER (from no.1)',/,
     $      ' ','             2) JT1 = JOINT NUMBER DEFINING ONE END OF 
     $THE MEMBER',/,
     $      ' ','             3) JT2 = JOINT NUMBER DEFINING THE OTHER E
     $ND OF MEMBER',/,
     $      ' ','             4) MI = MOMENT OF INERTIA, in4',/,
     $      ' ','             5) A = CROSS-SECTIONAL AREA, in2',/,
     $   ' ','             6) MP = PLASTIC MOMENT OF MEMBER, in-kips'/)
      DO 600 II=1,NM
         READ(IDSK,*) I,(MCON(I,J),J=1,2),SMA(I),AREA(I),PM(I)
         IF (STORE) WRITE(7,12) I,(MCON(I,J),J=1,2),SMA(I),AREA(I),PM(I)
  600 CONTINUE
   12 FORMAT(' ',3I5,3(F15.5,3X))
      L=0
      DO 700 I=1,JCT
      DO 700 J=1,3
  700 L=L+JTYPE(I,J)
      IF (IOPT .EQ. 1) WRITE(*,1010)
      IF (IOPT .EQ. 1) WRITE(*,'(A)')'      NUMBER OF JOINTS LOADED ?'
      READ(IDSK,*) LN
      IF (STORE) WRITE(7,'(I6)') LN 
      DO 720 I=1,L
  720 VL(I)=0.
      WRITE(8,305) JCT,NM,E
  305 FORMAT('%',5X,'* GENERAL DATA',/'%',10X,' NUMBER OF JOINTS',I16/
     ,'%',10X,' NUMBER OF MEMBERS',I15,/
     ,'%',10X,' MOD OF ELASTICITY',F15.1/'%'/'%')
      WRITE(8,505) 
  505 FORMAT('%',5X,' * DATA FOR JOINTS'/,
     ,'%',$10X,' JOINT   X-COORD   Y-COORD   DFX   DFY   DFZ',/)
      DO 506 I=1,JCT
  506 WRITE(8,507) I,(CORD(I,J),J=1,2),(JTYPE(I,J),J=1,3)
  507 FORMAT('%',I14,F12.2,F10.2,3I6)
      WRITE(8,605)
  605 FORMAT('%',/,'%',5X,' * DATA FOR MEMBERS'/,
     ,'%',$10X,' MEMBER   JT1     JT2       IXX      AREA        MP',/)
      DO 606 I=1,NM
  606 WRITE(8,607) I,(MCON(I,J),J=1,2),SMA(I),AREA(I),PM(I)
  607 FORMAT('%',I14,I9,I8,3F10.2)
      WRITE(8,705) 
  705 FORMAT('%',/'%',5X,' * DATA FOR LOADS'/,'%',
     $10X,' JOINT        PX        PY        PZ')
      IF (IOPT .EQ. 1) WRITE(*,1050)
 1050 FORMAT(//,'      ENTER "JT","FX","FY", and "FZ" in one line',/,
     $      ' ','                             (repeated LN times)',/,
     $      ' ','             1) JT = JOINT NUMBER',/,
     $      ' ','             2) FX = LOAD IN THE X-DIRECTION, kips',/,
     $      ' ','             3) FY = LOAD IN THE Y-DIRECTION, kips',/,
     $      ' ','             4) FZ = LOAD IN THE ROTATION DIRECTION, in
     $-kips',/)
      DO 770 I=1,LN
         READ(IDSK,*) JN,(OLEN(J),J=1,3)
         IF (STORE) WRITE(7,'(I6,5X,3F15.5)') JN,(OLEN(J),J=1,3)
         WRITE(8,706) JN,(OLEN(J),J=1,3)
  706    FORMAT('%',I14,F12.2,2F10.2)
         LL=0
         LJ=JN-1
         IF (LJ) 750,750,730
  730    DO 740 J=1,LJ
         DO 740 K=1,3
  740    LL=LL+JTYPE(J,K)
  750 DO 770 K=1,3
         IF (JTYPE(JN,K)) 770,770,760
  760    LL=LL+1
         VL(LL)=OLEN(K)
  770 CONTINUE
      WRITE(*,'(A)')' EXECUTION BEGINS.......'
      NCYCL=0
      CLG=0.
      DO 900 I=1,L
  900 CX(I)=0.
      WRITE(8,905) 0,0,0,0,(CX(J),J=1,3*JCT),(CM(K),K=1,2*NM),
     &(CT(I),I=1,NM)
      DO 270 I=1,NM
         J1=MCON(I,1)
         J2=MCON(I,2)
         X=CORD(J1,1)-CORD(J2,1)
         Y=CORD(J1,2)-CORD(J2,2)
         OLEN(I)=SQRT(X*X+Y*Y)
  270 CONTINUE
      DO 280 I=1,NM
         SF(2*I,2)=4.*E*SMA(I)/OLEN(I)
         SF(2*I-1,1)=SF(2*I,2)
         SF(2*I,1)=.5*SF(2*I,2)
         SF(2*I-1,2)=SF(2*I,1)
         SA(I)=E*AREA(I)/OLEN(I)
  280 CONTINUE
      M2=2*NM
      M3=3*NM
      DO 290 I=1,M2
  290 CM(I)=0.
      DO 291 I=1,NM
  291 CT(I)=0.
      NJ=0
      NK=0
      DO 295 I=1,L
      DO 295 J=1,M3
  295 A(I,J)=0.
      DO 450 J=1,JCT
         DO 440 M=1,NM
            NA=NJ
            IF (J-MCON(M,1)) 340,330,340
  330       JF=MCON(M,2)
            MJ=2*M-1
            MF=MJ+1 
            GOTO 360
  340       IF (J-MCON(M,2)) 440,350,440
  350       JF=MCON(M,1)
            MJ=2*M
            MF=MJ-1
  360       X=CORD(JF,1)-CORD(J,1)
            Y=CORD(JF,2)-CORD(J,2)
            D=SQRT(X*X+Y*Y)
            S=Y/D
            C=X/D
            NN=2*NM+M
            IF (JTYPE(J,1)) 380,380,370
  370       NA=NA+1
            A(NA,MJ)=S/D
            A(NA,MF)=A(NA,MJ)
            A(NA,NN)=-C
  380       IF (JTYPE(J,2)) 405,405,390
  390       NA=NA+1
            A(NA,MJ)=-C/D
            A(NA,MF)=A(NA,MJ)
            A(NA,NN)=-S
  405       IF (JTYPE(J,3)) 420,420,410
  410       NA=NA+1
            A(NA,MJ)=1.
  420       IF (NA-NK) 440,440,430
  430       NK=NA
  440    CONTINUE
         NJ=NK
  450 CONTINUE
  455 NCYCL=NCYCL+1
      DO 480 J=1,L
         DO 460 I=1,M2
            K=((I+1)/2)*2-1
            CSAT(I)=SF(I,1)*A(J,K)+SF(I,2)*A(J,K+1)
  460    CONTINUE
         DO 470 I=1,NM
            K=M2+I
            CSAT(K)=SA(I)*A(J,K)
  470    CONTINUE
      DO 480 I=1,L
         ASAT(I,J)=0.
      DO 480 K=1,M3
         ASAT(I,J)=ASAT(I,J)+A(I,K)*CSAT(K)
  480 CONTINUE
      DO 490 I=1,L
  490 ASAT(I,L+1)=VL(I)
      KJ=L+1
      DO 610 I=1,L
         IP1=I+1
         TEMP=ABS(ASAT(I,I))
         K=I
         DO 520 J=I,L
            IF (ABS(ASAT(J,I))-TEMP) 520,520,510
  510       K=J
            TEMP=ABS(ASAT(J,I))
  520    CONTINUE
         IF (K-I) 530,550,530
  530    DO 540 J=I,KJ
            TEMP=ASAT(I,J)
            ASAT(I,J)=ASAT(K,J)
            ASAT(K,J)=TEMP
  540    CONTINUE
  550    IF (ASAT(I,I)) 570,575,570
  570    TEMP=1./ASAT(I,I)
         DO 580 J=I,KJ
  580    ASAT(I,J)=ASAT(I,J)*TEMP
      DO 610 J=1,L
         IF (I-J) 590,610,590
  590    TEMP=ASAT(J,I)
         DO 608 K=IP1,KJ
  608    ASAT(J,K)=ASAT(J,K)-TEMP*ASAT(I,K)
  610 CONTINUE
      XLMT=1000.
      DO 620 I=1,L
         IF (ABS(ASAT(I,KJ))-XLMT) 620,620,625
  620 CONTINUE
      DO 650 I=1,M3
         CSAT(I)=0.
      DO 650 J=1,L
  650 CSAT(I)=CSAT(I)+A(J,I)*ASAT(J,L+1)
      DO 660 I=1,M2
         K=((I+1)/2)*2-1
         SATX(I)=SF(I,1)*CSAT(K)+SF(I,2)*CSAT(K+1)
  660 CONTINUE
      DO 670 I=1,NM
         K=M2+I
         SATX(K)=SA(I)*CSAT(K)
  670 CONTINUE
      DO 722 I=1,M2
         K=(I+1)/2
         ZERO=.001*PM(K)
         IF (ABS(SATX(I))-ZERO) 708,708,710
  708    ALG(I)=1.E10
         GOTO 722
  710    ALG(I)=(PM(K)-ABS(CM(I)))/ABS(SATX(I))
  722 CONTINUE
      SALG=1.E10
      DO 755 I=1,M2
         TEST=CM(I)*SATX(I)
         IF (TEST) 755,735,735
  735    IF (ALG(I)-SALG) 745,755,755
  745    SALG=ALG(I)
         NPH=I
  755 CONTINUE
      DO 765 I=1,M3
  765 SATX(I)=SALG*SATX(I)
      CLG=CLG+SALG
      DO 775 I=1,M2
  775 CM(I)=CM(I)+SATX(I)
      DO 780 I=1,NM
         K=M2+I
         CT(I)=CT(I)+SATX(K)
  780 CONTINUE
      DO 790 I=1,L
         ASAT(I,KJ)=ASAT(I,KJ)*SALG
         CX(I)=CX(I)+ASAT(I,KJ)
  790 CONTINUE
      I=(NPH+1)/2
      K=(NPH/2)*2-NPH
      IF (K) 792,793,793
  792 J=MCON(I,1)
      GOTO 799
  793 J=MCON(I,2)
  799 II = I 
      JJ = J
      WRITE(8,800) NCYCL,I,J,CLG
  800 FORMAT('%'/'%'/'%'/
     ,'%',5X,' * PLASTIC HINGE',I3,' FORMED IN MEMBER',I3
     ,' NEAR JOINT',I3,' WHEN LOAD FACTOR IS',F12.3/'%')
      CONTINUE
      WRITE(8,820)
  820 FORMAT('%',10X,' CUMULATIVE DEFORMATIONS')
      WRITE(8,821)
  821 FORMAT('%',15X,' JOINT    X-DISP       Y-DISP       ROTN')
      LL=0
      DO 830 I=1,JCT
         DO 827 J=1,3
            IF (JTYPE(I,J)) 825,825,826
  825       CSAT(J)=0.
            GOTO 827
  826       LL=LL+1
            CSAT(J)=CX(LL)
  827    CONTINUE
         WRITE(8,840) I,(CSAT(J),J=1,3)
  830 CONTINUE
  840 FORMAT('%',I19,3F13.5)
      WRITE(8,850)
  850 FORMAT('%'/,'%',10X,' CUMULATIVE MOMENTS')
      WRITE(8,851)
  851 FORMAT('%',15X,' MEMBER       END MOMENTS       ',
     $'    JOINTS     PLASTIC MOM')
      DO 860 I=1,NM
         K=2*I-1
         WRITE(8,870) I,CM(K),CM(K+1),MCON(I,1),MCON(I,2),PM(I)
  860 CONTINUE
  870 FORMAT('%',I19,F14.2,F11.2,I7,' AND',I2,F14.2)
      WRITE(8,880)
  880 FORMAT('%',/,'%',10X,' CUMULATIVE TENSION FORCES')
      WRITE(8,885)
  885 FORMAT('%',15X,' MEMBER     TENSION')
      DO 890 I=1,NM
         K=M2+I
         WRITE(8,901) I,CT(I)
  890 CONTINUE
  901 FORMAT('%',I19,F15.2)
      WRITE(8,905) NCYCL,II,JJ,CLG,(CX(J),J=1,3*JCT),(CM(K),K=1,2*NM),
     &(CT(I),I=1,NM)
  905 FORMAT('%',/I4,I4,I4,1000E16.4)
      ITEST=((NPH/2)*2)-NPH
      IF (ITEST) 910,920,920
  910 SF(NPH+1,2)=.75*SF(NPH+1,2)
      SF(NPH+1,1)=0.
      SF(NPH,1)=0.
      SF(NPH,2)=0.
      GOTO 455
  920 SF(NPH-1,1)=.75*SF(NPH-1,1)
      SF(NPH-1,2)=0.
      SF(NPH,1)=0.
      SF(NPH,2)=0.
      GOTO 455
  575 WRITE(8,576)
  576 FORMAT('%',5X,' * DIVISION BY ZERO IN SOLUTION OF EQUATION'/,'%')
      GOTO 940
  625 WRITE(8,923) XLMT,NCYCL
  923 FORMAT('%',/,'%',5X,' *** DEFORMATIONS LARGER THAN ',F8.1,
     $' IN CYCLE NO',I4/,'%')
  940 WRITE(8,950) JFN
  950 FORMAT('%',/,'%',5X,' ANALYSIS COMPLETED FOR FRAME NO',I3/)
      WRITE(*,'(A)')'              ............ EXECUTION TERMINATED'
      STOP 
      END
