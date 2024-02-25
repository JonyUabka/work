*
*+------------------------------------------------------------
*+    程序功能：计算滤波因子
*+    参数：
*+        LMS--     滤波算子长度（ms）
*+        ISI--     采样间隔（ms）
*+        JOB--     滤波因子类型
*+                  =1  低通滤波
*+                  =2  高通滤波
*+                  =3  陷频滤波
*+                  =4  带通滤波
*+        RIT1--    带通滤波参数1
*+        RIT2--    带通滤波参数2
*+        RIT3--    带通滤波参数3
*+        RIT4--    带通滤波参数4
*+        FF--      存放滤波算子
*+------------------------------------------------------------
*
* 
      SUBROUTINE CALFIL(LMS,ISI,JOB,RIT1,RIT2,RIT3,RIT4,FF)
      REAL  FF(*)
* 
      NO=LMS/(2*ISI) 
      NT=2*NO+1 
      TAU=ISI 
      JEC=0 
* 
      IF ( JOB.NE.3 ) GOTO 400 
* 
*+    JOB = 3, computing Notch operator NH(f1,f2,f3). 
*+    Adjusting the originally defined frequencies (f1,f2,f3) 
*+    to the form of (f1,f2,f2,f3). 
* 
      IF ( RIT4.NE.0 ) GOTO 211 
         RIT4=RIT3 
         RIT3=RIT2 
         GOTO 211 
* 
  400 IF ( JOB.EQ.4 ) GOTO 210 
* 
*+    For High and Low Pass operator,the frequency is (0,0,f3,f4). 
* 
         RIT3=RIT1 
         RIT4=RIT2 
         RIT1=0 
         RIT2=0 
* 
      IF ( JOB.NE.2 ) GOTO 210 
* 
  211 JEC=1 
* 
*+    JEC -- a flag used to compute high pass and notch operator 
* 
  210 CONTINUE 
* 
*+    For Band Pass operator,the form of frequency is (f1,f2,f3,f4) 
*+ 
*+    TMI  = ( f4 - f3 ) / 2.0 
*+    TMA  = ( f4 + f3 ) / 2.0 
*+    TMII = ( f2 - f1 ) / 2.0 
*+    TMAA = ( f2 + f1 ) / 2.0 
*+ 
      TMI=(RIT4-RIT3)/2.0 
      TMA=(RIT3+RIT4)/2.0 
      TMII=(RIT2-RIT1)/2.0 
      TMAA=(RIT1+RIT2)/2.0 
      FIM=0.0 
      PI=3.141593 
      K=NO+2 
      IC=1 
      TI=TMI 
      TA=TMA 
CLIP2004.10.06      TER=2.0*PI*TAU*0.001 
C      TER=2.0*PI*TAU*0.001/1000.
      TER=2.0*PI*TAU*0.001
  311 CONTINUE 
      B=TER*TA 
      BETA=TER*TI 
* 
*+    B    = PI * TAU * ( f4 + f3 )  or   PI * TAU * ( f2 + f1 ) 
*+    BETA = PI * TAU * ( f4 - f3 )  or   PI * TAU * ( f2 - f1 ) 
* 
      DO 10 I=1,NO 
* 
         Z=I 
         X=B*Z 
         Y=BETA*Z/2.0 
         SINB=SIN(X) 
         SINC=SIN(Y) 
* 
*+       SINB = sine [ PI*TAU*i*(f4+f3) ] 
*+       SINC = sine [ PI*TAU*i*(f4-f3)/2.0 ] 
* 
         SIN2=SINC*SINC 
         COE=BETA*BETA*Z*Z*Z/4.0*PI 
* 
*+       COE is the coefficient 
* 
         FF(K)=(SINB*SIN2*1.0)/COE 
         K=K+1 
* 
   10 CONTINUE 
* 
      GOTO (320,321),IC 
* 
*+    The symmetrical point of the operator is TAU * ( f4 + f3 ) 
*+    If f1 or f2 is not zero,computing again. 
* 
  320 FF(NO+1)=(B*1.0)/PI 
      IF (TMAA.LE.0.0) GOTO 330 
      IC=2 
      TI=TMII 
      TA=TMAA 
      K=1 
      GOTO 311 
* 
*+    Band Pass operator 
* 
  321 CONTINUE 
      DO 501 L=1,NO 
         JJ=NO+L+1 
         FF(JJ)=FF(JJ)-FF(L) 
  501 CONTINUE 
      FIM=(B*1.0)/PI 
* 
*+    Completing the operator by sending one side of the operator 
*+    to the other side according to the symmetrical point. 
* 
  330 CONTINUE 
      DO 502 L=1,NO 
         II=NO+L+1 
         JJ=NO-L+1 
         FF(JJ)=FF(II) 
  502 CONTINUE 
* 
      FF(NO+1)=FF(NO+1)-FIM 
      IF ( JEC.EQ.0 ) RETURN 
* 
*+    Substraction of a pulse to compute High-Pass and Notch operator 
* 
      FF(NO+1)=FF(NO+1)-1.0 
      DO 506 J=1,NT 
         FF(J)= -FF(J) 
  506 CONTINUE 
* 
*      write(*,*)"FF==",(FF(IJ),IJ=1,NT)
      RETURN 
      END 

*
*+------------------------------------------------------------
*+    程序功能：对地震道进行滤波
*+    参数：
*+        filt --  滤波因子
*+        nt   --  滤波因子长度
*+        LSAMP--  地震道长
*+        trace--  地震道
*+------------------------------------------------------------
*

      subroutine filter_exec( FILT,NT,LSAMP, TRACE)
      real TRACE(*), FILT(*), WKPADI(LSAMP)
C

*
*+    将输入的地震道trace传送道工作缓冲区
*
*      write(*,*)"enter filter LSAMP==",LSAMP
      DO IK=1,LSAMP
         WKPADI(IK)=TRACE(IK)
      ENDDO

*         write(*,*) "trace* = ", (trace(IK),IK=10,20)
*         write(*,*) "WKPADI* = ", (WKPADI(IK),IK=10,20)

*         write(*,*) "FILT* = ", (FILT(IK),IK=1,NT)
         IFLAG=-1
         ISHIFT=NT/2+1
*
*+    褶积
*
*         write(*,*)"LSAMP==",LSAMP,"NT==",NT
*         write(*,*)"FILT==",(FI(IJ),IJ=1,NT)
         CALL CONVOX(TRACE,LSAMP,WKPADI,LSAMP,FILT,NT,
     +IFLAG,ISHIFT)
*      write(*,*) "trace = ", trace(101)
       RETURN
       END

*
*+------------------------------------------------------------
*+    程序功能：相关或褶积运算
*+    参数：
*+        Y--      输出数组
*+        LY--     输出数组长度
*+        X--      输入数组
*+        LX--     输入数组长度
*+        F--      滤波算子
*+        LF--     滤波算子长度
*+        IFLAG--  相关或褶积
*+                 >=0，相关
*+                 其它，褶积
*+        ISHIFT-- 時移量
*+------------------------------------------------------------
*

      SUBROUTINE CONVOX(Y,LY,X,LX,F,LF,IFLAG,ISHIFT)
      REAL  Y(*),X(*),F(*)
* 
*     Firstly calculating the shifting samples 
* 
*      write(*,*) "X==",(X(IK),IK=1,LX,50)
*      write(*,*) "LY,LX,LF==",LY,LX,LF,IFLAG,ISHIFT
*      write(*,*) "X==",(X(IK),IK=1,LX)
*      write(*,*) "F==",(F(IK),IK=1,LF)
*       return
      IF (IFLAG.GE.0) THEN 
*     -------------------- 
         II=ISHIFT 
      ELSE 
*     ---- 
         II=ISHIFT-LF+1 
         DO 5 I=1,LF/2 
            FC=F(I) 
            F(I)=F(LF-I+1) 
            F(LF-I+1)=FC 
    5    CONTINUE 
      ENDIF 
*     ----- 
*     Convoluting 
*      write(*,*) "F==",(F(IK),IK=1,LF)

      DO  100  I = 1 , LY 
* 
          ASUM=0.0
          IS=1 
          IFL=LF 
          IA=II 
* 
*     Here IA is start sample of X , IS and IFL are start and end 
*     samples of F respectively. 
* 
          IF (IA.LT.0) THEN 
*         ----------------- 
             IS=-IA+1 
             IA=0 
          ELSE 
*         ---- 
             IF (IA+LF.GT.LX) THEN 
*            --------------------- 
                IFL=LX-IA 
             ENDIF 
*            ----- 
          ENDIF 
*         ----- 
*         write(*,*) "IS,IFL,IA=",IS,IFL,IA
*         return
              IF (IS.GT.IFL) STOP
*
*        write(*,*) "X==",(X(IK),IK=1,LX,50)
          DO 50 J=IS,IFL
             IA=IA+1 
*             write(*,*)"X(IA),F(J)=",X(IA),F(J),IA,J
             ASUM=ASUM+X(IA)*F(J)
   50     CONTINUE 
          II=II+1 
          Y(I)=ASUM
  100 CONTINUE 
*      WRITE(*,*) "Y==",(Y(IK),IK=1,LY,50)
      RETURN 
      END
