*
* @file prdecon.f
* @brief 完成预测反褶积功能
*
*
*
* @brief 完成预测反褶积模块的参数分析，工作缓冲区长度计算等功能。
* @param[in] IPAR 界面参数
* @param[in] DT  采样间隔
* @param[in] NT  道长
* @param[out] LC  工作缓冲区长度
*
      SUBROUTINE PRDECON_AM (address1,params,IPAR, DT, NT, LC )
      integer params(12)
      integer*8 address1

      DIMENSION  IPAR(*)
      DIMENSION  IPAM(800),IPON(30)
      INTEGER DT
      COMMON /PRDECON_PAR1/IPAM
      COMMON /PRDECON_PAR2/IPON
*

c      write(*,*) "DT=",DT,"NT=",NT
      I1=23            
      L1=800          
      LR1=0                                                      
      LR2=0
      I2=I1+L1
c c     write(*,*) "I2==",I2
*
* @brief   对界面参数进行处理，形成公共参数表。
* @param[in] IPAR 界面参数
* @param[out] IPAM 公共参数表
* @param[in] DT  采样间隔
* @param[in] NT  道长
*
      CALL PRDECON_PARM ( address1,params,IPAR,IPAM,DT,NT )
*
      L2=LR1       
      I3=I2+L2     
      L3=LR2                                                     
*
      I4=I3+L3
      L4=LR2
*
      K1=2
      LK1=128
      K2=K1+LK1  
      LK2=NT/DT
c      write(*,*) "LK2==",LK2
      K3=K2+LK2 
      K4=K3+LK2
      K5=K4+LK2
      K6=K5+257 
      K7=K6+LK2 
*
C      write(*,*) "L1,L2,L3,L4==",L1,L2,L3,L4
C      write(*,*) "LK2,LK2==",LK2,LK2
      LP=23+L1+L2+L3+L4
C      write(*,*) "LP==",LP
*
* @brief  计算工作缓冲区长度。
*
      LC=1+LK1+5*LK2+257


*
* @brief  将工作缓冲区起始地址存放在IPON数组中。
*
c      write(*,*) "LC1==",LC

      IPON(3)= I1
      IPON(4)= I2
      IPON(5)= I3
      IPON(6)= I4
c      write(*,*) "K1,K2,K3,K4,K5,K6,K7=",K1,K2,K3,K4,K5,K6,K7
      IPON(7)= K1
      IPON(8)= K2
      IPON(9)= K3
      IPON(10)= K4
      IPON(11)=K5
      IPON(12)=K6
      IPON(13)=K7
C      write(*,*) "out of decon!!!"

*
*
      RETURN
      END
*
*
*
*
* @brief 完成预测反褶积模块的核心算法。
* @param[in] KBUF 工作缓冲区
* @param[in] LIBM 切除数据表
* @param[in/out] IH  道头
* @param[in/out] TRACE  地震道
*
      SUBROUTINE PRDECON_PM (KBUF, LIBM, IH, TRACE )
      DIMENSION  KBUF(*),IH(*),TRACE(*)
      DIMENSION  IPAM(800),LIBM(*),IPON(30)
C      character cmname*80
      COMMON /PRDECON_PAR1/IPAM
      COMMON /PRDECON_PAR2/IPON
*

*+      Memory distribution table
*+     +----------------------------------------------+
*+     | I/O registers                       len.=2   |
*+     +----------------------------------------------+
*+     | IPAM                                len.=200 |
*+     |     private buffer for decoding parameters   |
*+     +----------------------------------------------+
*+     | LIBM                                len.=L2  |
*+     |     private buffer for mute library          |
*+     +----------------------------------------------+
*=     | IM                                  len.=L3  |
*+     |     private buffer for mute function         |
*+     +----------------------------------------------+
*+     | IX                                  len.=L4  |
*+     |     private buffer for mute function         |
*+     +----------------------------------------------+
*+     +----------------------------------------------+
*+     | IH                                  len.=128 |
*+     |     common buffer for trace header           |
*+     +----------------------------------------------+
*+     | E   common buffer for trace data             |
*+     +----------------------------------------------+
*+     | the following are common working areas       |
*+     |     ZSC                                      |
*+     |     ZS                                       |
*+     |     IV                                       |
*+     |     COEF                                     |
*+     |     TRF                                      |
*+     +----------------------------------------------+
*
      I1=IPON(3)
      I2=IPON(4)
      I3=IPON(5)
      I4=IPON(6)
      K1=IPON(7)
      K2=IPON(8)
      K3=IPON(9)
      K4=IPON(10)
      K5=IPON(11)
      K6=IPON(12)
      K7=IPON(13)
*
      LTR=IH(10)/IH(9)
*      MAX=IABS(IH(14))
*
      IF (IH(11).EQ.2) THEN
*     -----------------------------------
c       write(*,*) "LTR===",LTR
       DO 15 I=1,LTR
         TRACE(I)=0
   15  CONTINUE
       RETURN
      ENDIF
*     -----
*
       IF (IPAM(21).NE.0) THEN
*
* @brief 计算当前道的切除时间。
* @param[in] IPAM 公共参数表
* @param[in] LIBM 存放切除数据
* @param[in] IH  道头
*
        CALL PRDECON_MUTE ( IPAM ,LIBM ,IH )
c        write(*,*)"IH(20),IPAM(83),IPAM(86)=",IH(20),IPAM(83),IPAM(86)
        DO IK=1,IPAM(83)/IH(9)
           TRACE(IK)=0
        ENDDO
*
       ELSE
*     ----
       IPAM(83)=IH(19)
       IPAM(86)=IH(7)
       ENDIF
*     -------
* @brief 算法核心程序，完成预测算子计算及反褶积功能。
* @param[in] IPAM 公共参数表
* @param[in/out] TRACE  地震道
* @param[] KBUF(K3)  工作缓冲区，存放中间结果
* @param[] KBUF(K4)  工作缓冲区，存放计算参数
* @param[] KBUF(K6)  工作缓冲区，存放中间结果
* @param[] KBUF(K7)  工作缓冲区，存放输入地震道
* @param[in/out] IH  道头
**
C      write(*,*)"K3,K4,K5,K6,K7=",K3,K4,K5,K6,K7
      CALL PRDECON_CORE (IPAM,TRACE,KBUF(K3),KBUF(K4),
     +               KBUF(K6),KBUF(K7),IH )
C      write(*,*)"out of PRDECON_CORE!!!"
*
* @brief 根据选件是否做振幅调整。
* @param[in] IPAM 公共参数表
* @param[in/out] TRACE  地震道
* @param[] KBUF(K3)   工作缓冲区
* @param[] KBUF(K4)   工作缓冲区
* @param[] KBUF(K5)   工作缓冲区
* @param[in/out] IH   道头
*
      CALL PRDECON_ADJUST (IPAM,TRACE,KBUF(K3),
     +             KBUF(K4),KBUF(K5),IH)

      RETURN
      END
*
*
* @brief   对界面参数进行处理，形成公共参数表。
* @param[in] IPAR 界面参数
* @param[out] IPAM 公共参数表
* @param[in] DT  采样间隔
* @param[in] NT  道长
*
      SUBROUTINE PRDECON_PARM (address1,TMTABLE,IPAR,IPAM,DT,NT)
      DIMENSION IPAR(*),IPAM(*)
      INTEGER DT ,SPACEFLAG
      integer  TMTABLE(12)
      integer*8 address1
      CHARACTER*150 bbb
      integer intnum(150)
      integer  temptable1(4),temptable2(4)
      integer  temptable3(4),temptable4(4)


*
*=========================================================*
*    Parameter Table in IPAM                              *
*                                                         *
*    IPAM (01-05) : W1 the first window limits            *
*         (06-10) : W2 the second window limits           *
*         (11-15) : PD the period lengths                 *
*         (16-20) : CF1 the coefficient                   *
*         (   21) : i, number of mute library             *
*         (   22) : =1 (WB)                               *
*         (   23) : TI limit of the shallow zone          *
*         (   24) : LW                                    *
*         (   25) : DEB                                   *
*         (   26) : L/2                                   *
*         (   27) : =1( A )  or =0( B )                   *
*         (28-31) : frey                                  *
*         (   32) : =1(LIST)                              *
*         (   33) : frey                                  *
*         (   36) : =1 (if PD)                            *
*         (   37) : =1 (if TI)                            *
*         (   38) : L                                     *
*         (   40) : number of window                      *
*                                                         *
*         (201)   HD1                                     *
*                                                         *
*                                                         *
*                                                         *
*=========================================================*


C      call GetParam_int(params,"spaceflag"//CHAR(0),SpaceFlag)
C      write(*,*)"SpaceFlag==",SpaceFlag
c      write(*,*) "TMTABLE=",(TMTABLE(ik),ik=1,12)
      SpaceFlag = IPAR(5)
      if(SpaceFlag.ne.1) then
!         ColNUM=12
!         Icon=1
!         DO IK=1,ColNUM
!           call get_table_param_columnvalues_int(params,
!     &     "veltable1"//CHAR(0),tabname1(IK),intnum(icon),intcount)
!           icon=icon+intcount
!         ENDDO
!         write(*,*)"1intnum ,intcount==",intnum ,intcount


CCCCCCCCCCCCCCCCCC 郭利荣修改 获取不空变的时窗 2012.4.19 CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

           do ik=1,4
              ikk=(ik-1)*3
              intnum(ik)=TMTABLE(1+ikk)
              intnum(4+ik)=TMTABLE(2+ikk)
              intnum(8+ik)=TMTABLE(3+ikk)
           enddo
c           write(*,*) "intnum=",(intnum(ik),ik=1,12)




c  为什么要将intcount赋为1呢？
           intcount =1
!       write(*,*)"1intnum ,intcount==",intnum ,intcount
!       STOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



         IWINCOUNT=0
         IERR=0
         DO IK=3,1,-1
            IF(intnum(IK*4*intcount).NE.0) THEN
               IWINCOUNT=IWINCOUNT+1
               DO IJ=1,3*intcount-1
                  IF(intnum(IK*4*intcount-IJ).EQ.0) THEN

                     bbb="The gate parameter error!"
*                     call WriteLog(ADDRESS1, 0,bbb)
                     IERR=IERR+1
                  ENDIF
               ENDDO
            ELSE
               DO IJ=1,4*intcount-1
                  IF(intnum(IK*4*intcount-IJ).NE.0) THEN

                     bbb="The gate parameter error!"
*                     call WriteLog(ADDRESS1, 0,bbb)
                     IERR=IERR+1
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         IF(IWINCOUNT.EQ.0) THEN
            bbb="The gate parameter error!"
*            call WriteLog(ADDRESS1, 0,bbb)
            STOP
         ENDIF
         DO IK=1,IWINCOUNT
            imm=(IK-1)*4*intcount+1
            IF(intnum(imm+intcount).LE.intnum(imm)) THEN
            bbb="The gate parameter reversed!"
*            call WriteLog(ADDRESS1, 0,bbb)
               IERR=IERR+1
            ENDIF
         ENDDO
         IF(IERR.NE.0) STOP
C
      else
         ColNUM=14
         Icon=1
         DO IK=1,ColNUM
c           call get_table_param_columnvalues_int(params,
c     &     "veltable"//CHAR(0),tabname2(IK),intnum(icon) ,intcount)
           icon=icon+intcount
         ENDDO
!         write(*,*)"2intnum ,intcount==",intnum ,intcount
         IWINCOUNT=0
         IERR=0
         DO IK=3,1,-1
            IF(intnum(IK*4*intcount+2*intcount).NE.0) THEN
               IWINCOUNT=IWINCOUNT+1
               DO IJ=1,3*intcount-1
                  IF(intnum(IK*4*intcount+2*intcount-IJ).EQ.0) THEN
                     bbb="The gate parameter error!"
*                     call WriteLog(ADDRESS1, 0,bbb)
                     IERR=IERR+1
                  ENDIF
               ENDDO
            ELSE
               DO IJ=1,4*intcount-1
                  IF(intnum(IK*4*intcount+2*intcount-IJ).NE.0) THEN
                     bbb="The gate parameter error!"
*                     call WriteLog(ADDRESS1, 0,bbb)
                     IERR=IERR+1
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         IF(IWINCOUNT.EQ.0) THEN
            bbb="The gate parameter error!"
*            call WriteLog(ADDRESS1, 0,bbb)
            STOP
         ENDIF
         DO IK=1,IWINCOUNT
            DO IJ=1,intcount
               imm=(IK-1)*4*intcount+2*intcount+IJ
               IF(intnum(imm+intcount).LE.intnum(imm)) THEN
                  bbb="The gate parameter reversed!"
*                  call WriteLog(ADDRESS1, 0,bbb)
                  IERR=IERR+1
               ENDIF
            ENDDO
         ENDDO
         IF(IERR.NE.0) STOP
      endif


C      call GetParam_int(params,"length"//CHAR(0),ipar(2))
C      call GetParam_int(params,"muteflag"//CHAR(0),muteflag)
C      write(*,*) "muteflag==",muteflag
C      if(MuteFlag.eq.1) then

*       string MuteFile = param->GetParameterItem("mutefile")->ValueToString();
C       endif
C      call GetParam_int(params,"adjustflag"//CHAR(0),ipar(1))

C      if(ipar(1).eq.1) then
C          call GetParam_int(params,"lsmoth"//CHAR(0),ipar(3))
C          call GetParam_int(params,"shift"//CHAR(0),ipar(4))
C      endif
c      write(*,*) "ipar==",(ipar(ij),ij=1,6)
c      write(*,*) "read parameter end!!!!!!!!!!!!!!"

*
c      write(*,*) "enter decon1=="
c      write(*,*)"2window number =====",IWINCOUNT
c      write(*,*)"2location number =====",intcount
      DO 157 I=1,800
        IPAM(I)=0
  157 CONTINUE
      IPAM(44)=40
      IPAM(24)=10
      ISI=DT
      NE1=NT/DT
      IPAM(27)=1
        IPAM(28)=250/DT
        LENMUT = 0
      LTAB=20
      LWFLAG=0
      I=0
      IERR=0
      IREF=0
      IPAM(201)=5
      IPAM(202)=intcount
      IPAM(203)=IWINCOUNT
      IPAM(301)=1
      IPAM(302)=10000
      IFUNC=0
      IWIN=0

      IPAM(21)=IPAR(6)
      IPAM(23)=1
      IPAM(24)=IPAR(3)/(ISI*2)
      IPAM(25)=IPAR(4)
      IPAM(26)=IPAR(2)/(ISI*2)
      IPAM(27)=IPAR(1)
      IPAM(36)=1
      IPAM(38)=IPAM(26)*2+1
      IPAM(40)=IWINCOUNT
C      IREF=IPAR(7)/ISI+IPAM(26)
C
      ICON=3
      IPAN=4
      IDIS=5
      NICO=2+IPAN*ICON
      IREF=0
C      IADDRESS=0
      DO IJ =1,intcount
         IADDRESS=-intcount
         IF(SpaceFlag.EQ.1) THEN
         IADDRESS=IADDRESS+intcount
         IPAM(301+25*(IJ-1))=intnum(IJ)
         IADDRESS=IADDRESS+intcount
         IPAM(302+25*(IJ-1))=intnum(IJ+IADDRESS)
         ENDIF
         DO IK =1,IWINCOUNT
            intnum(IJ+IADDRESS+3*intcount)=
     &            intnum(IJ+IADDRESS+3*intcount)/ISI+IPAM(26)
            IF(IREF.LT.intnum(IJ+IADDRESS+3*intcount))
     &         IREF=intnum(IJ+IADDRESS+3*intcount)
            DO IM=1,IPAN
               IADDRESS=IADDRESS+intcount
               IPAM(305+25*(IJ-1)+IK+(IM-1)*IDIS)= intnum(IJ+IADDRESS)
            ENDDO
         ENDDO
      ENDDO
C      write(*,*) "IREF==",IREF
C      DO II=1,60
C      write(*,*) "II=",II,(IPAM(II*10-10+IJ),IJ=1,10)
C      ENDDO

C      IREF=IREF/ISI+IPAM(26)
C
      IQMAX=0
      IF(IPAM(36).EQ.1) IQMAX=IREF
      IPAM(41)=IQMAX+IPAM(26)+1
C      IPAM(161)=0
      J=IPAM(40)
      DO 102 IJ=1,J
         IPAM(160+J)=IPAM(41)
102   CONTINUE

c      write(*,*) "out of decon1=="
C      LIBM(1)=2
C      LIBM(2)=2
C      LIBM(3)=100
C      LIBM(5)=0
C      LIBM(6)=1001
C      LIBM(7)=1005
C      LIBM(8)=1002
C      LIBM(9)=1006
C      LIBM(10)=20
C      LIBM(11)=24
C      LIBM(12)=2
C      LIBM(13)=3
C      LIBM(20)=500
C      LIBM(21)=2000
C      LIBM(22)=50
C      LIBM(23)=4000
C      LIBM(24)=400
C      LIBM(25)=1000
C      LIBM(26)=3000
C      LIBM(27)=40
C      LIBM(28)=2000
C      LIBM(29)=5000

      return
      END

*
* @brief 计算当前道的切除时间。
* @param[in] IPAM 公共参数表
* @param[in] LIBM 存放切除数据
* @param[in] IH  道头
*
      SUBROUTINE PRDECON_MUTE (IPAM,LIBM,IH)
      DIMENSION  IPAM(*),LIBM(*),IH(*)
      DIMENSION IBINF(50),IBSUP(50),NL(50),IP(50)
      INTEGER   TMU,TMU1,TMU2,KTAP
*
*
      IF ( IPAM(22).EQ.1 ) THEN
*     -------------------------
      ICDP=IH(19)
      IQ=IH(19)
      ELSE
*     ----
       IQ=0
       KVAL=LIBM(1)           !空变道头字
      ICDP=IH(KVAL)
      ENDIF
*     -----
      IXX=IH(20)
      NU=LIBM(2)              !空变点数
      KTAP=LIBM(3)            !切除斜坡
      IFLAG=LIBM(5)           !标志
C      write(*,*) "NU,KTAP,IFLAG,KVAL=",NU,KTAP,IFLAG,KVAL
*
       DO 1 I=1,NU
       IBINF(I)=LIBM(5+I)     !每组起始CDP
       IBSUP(I)=LIBM(5+I+NU)  !每组终止CDP
       IP(I)=LIBM(5+2*NU+I)   !每组切除值存放起始地址
       NL(I)=LIBM(5+3*NU+I)   !每组TMU-TX对数
c       write(*,*) "I,IBINF(I),IBSUP(I)=",I,IBINF(I),IBSUP(I)
c       write(*,*) "I,IP(I),NL(I)=",I,IP(I),NL(I)
  1    CONTINUE
*
*
       I=1
  25   CONTINUE
*
c       write(*,*) "ICDP=",ICDP
       IF(ICDP.GT.IBSUP(I).AND.ICDP.LT.IBINF(I+1)) THEN
       IPR1=IP(I)
       IPR2=IP(I+1)
       NBC1=NL(I)
       NBC2=NL(I+1)
*
*
       K1=IPR1
       K2=IPR1+NBC1
       K3=IPR2
       K4=IPR2+NBC2
!       write(*,*) "K1,K2,K3,K4=",K1,K2,K3,K4
!       write(*,*) "NBC1,NBC2,IXX=",NBC1,NBC2,IXX
*
*
      CALL PRDECON_MUTE_INT(LIBM(K1),LIBM(K2),NBC1,IXX,TMU1)
      CALL PRDECON_MUTE_INT(LIBM(K3),LIBM(K4),NBC2,IXX,TMU2)
*
         IF(IBSUP(I).EQ.0) THEN
         VAL1=ICDP-IBINF(I)
         VAL2=IBINF(I+1)-IBINF(I)
         ELSE
         VAL1=ICDP-IBSUP(I)
         VAL2=IBINF(I+1)-IBSUP(I)
         ENDIF
          TMU=TMU1+(TMU2-TMU1)*VAL1/VAL2
*
* @brief  IPAM(83)  mute value
* @brief  IPAM(86)  mute taper
*
          IPAM(83)=TMU
          IPAM(86)=KTAP
        IF ( IFLAG.NE.0.AND.IPAM(22).NE.1)THEN
        IPAM(83)=IPAM(83)+IQ
        ENDIF
          RETURN
      ENDIF
*
        IF(ICDP.LE.IBSUP(I)) THEN
*
          IPR1=IP(I)
          NBC1=NL(I)
          K1=IPR1
          K2=IPR1+NBC1
c       write(*,*) "2K1,K2=",K1,K2
c       write(*,*) "NBC1,IXX=",NBC1,IXX
          CALL PRDECON_MUTE_INT(LIBM(K1),LIBM(K2),NBC1,IXX,TMU)
*
c          write(*,*) "TMU=",TMU,"KTAP=",KTAP,"IQ=",IQ
          IPAM(83)=TMU
          IPAM(86)=KTAP
      IF(IFLAG.NE.0.AND.IPAM(22).NE.1) THEN
       IPAM(83)=IPAM(83)+IQ
      ENDIF
          RETURN
        ENDIF
*
        I=I+1
c        IF( IBINF(I).EQ.0) GO TO 170
        GO TO 25
*
*
c 170   write(*,*) 'Current CMP number is out of limit defined'
c     1 //'in MU data table.'

        END
*     -----
*
*
*
* @brief 对切除库进行处理，获得空变点的切除时间。
* @param[in] M 切除时间数组
* @param[in] IX 偏移距数组
* @param[in] N 时间-偏移距对
* @param[in] ID 当前道偏移距
* @param[in] ITM  插值出的当前偏移距的切除时间
*
      SUBROUTINE PRDECON_MUTE_INT (M,IX,N,ID,ITM)
      DIMENSION  M(*),IX(*)
c      write(*,*) "ID=",ID,"N==",N
c      write(*,*) "M=",(M(IK),IK=1,N)
c      write(*,*) "IX=",(IX(IK),IK=1,N)
      IF (ID.LE.IX(1)) THEN
         ITM=M(1)
         RETURN
      ENDIF
      DO 5 I=2,N
         IF (ID.LT.IX(I)) GOTO 10
    5 CONTINUE
      ITM=M(N)
      RETURN
   10 VAL1=M(I)-M(I-1)
      VAL2=IX(I)-IX(I-1)
      ITM=M(I-1)+(ID-IX(I-1))*VAL1/VAL2
      RETURN
      END
*
*
* @brief 算法核心程序，完成预测算子计算及反褶积功能。
* @param[in] IPAM 公共参数表
* @param[in/out] E  地震道
* @param[] ZS  工作缓冲区，存放中间结果
* @param[] IV  工作缓冲区，存放计算参数
* @param[] ZSC  工作缓冲区，存放中间结果
* @param[] TRF  工作缓冲区，存放输入地震道
* @param[in/out] IH  道头
*
      SUBROUTINE PRDECON_CORE(IPAM,E,ZS,IV,ZSC,TRF,IH)
      DIMENSION  IPAM(*),E(*),ZS(*),IV(*)
      DIMENSION  LOP(3),IDEL(3),INDEX1(10),T11(10),T21(10)
      DIMENSION  TRF(*),IH(*),ZSC(*),IQZ(5),IGZ(5),MMZ(5)
      DIMENSION  TT1(10),TT2(10),TT3(10),TT4(10)
      INTEGER    TT1,TT2,TT3,TT4,TM,TMU,TI,C,B,T11,T21,TM1
       INTEGER    TI1
*
*
*+   calculate window
*
      IHD=IPAM(201)
      IFUNC=IPAM(202)
      IWIN=IPAM(203)
      INHED=IH(IHD)
c      write(*,*)"INHED==",INHED
*
* @brief 通过插值计算算子参数。
* @param[in] INHED  当前空变点值
* @param[in] IFUNC  空变组数
* @param[in] IWIN  时窗个数
* @param[in/out] IPAM  公共参数表
* @param[in] IPAM（301） 空变参数
*
      CALL PRDECON_CORE_CALPARM(INHED,IFUNC,IWIN,IPAM,IPAM(301))

C      print *,"IPAM(6,22,86)=",IPAM(6),IPAM(22),IPAM(86)
C      print *,"IPAM(81,83,40)=",IPAM(81),IPAM(83),IPAM(40)
      IE=IH(16)
      KTAP=IPAM(86)
      NLE=IH(1)/IH(9)
      NE1=IH(10)/IH(9)
      TM=IPAM(81)/IH(9)*IH(9)
      K=IH(1)-IH(9)*11
*     IF(TM.GT.K)GO TO 900
      TMU=IPAM(83)-IE
      TI=TMU/IH(9)*IH(9)
      LENGTH=IH(1)/IH(9)
C      write(*,*) "TMU,TI,LENGTH=",TMU,TI,LENGTH
      N=IPAM(40)
      ISI=IH(9)
       IQQ=0
       IF( IPAM(22).EQ.1) THEN
      IQQ=IH(19)/ISI
       ENDIF
      NPM=IH(4)
*
      DO 10 I=1,5
         TT1(I)=IPAM(I)/ISI
         TT2(I)=IPAM(I+5)/ISI
         IQZ(I)=IPAM(I+10)
  10  CONTINUE
*
C      print *,"IPAM(36,26)=",IPAM(36),IPAM(26)
      IF ( IPAM(36).NE.1 ) THEN
*       -------------------------
       DO 11 I=1,5
      IQZ(I)=IQQ
 11     CONTINUE
      ENDIF
*       -----
*
C      print *,"IQZ(1-3)=",IQZ(1),IQZ(2),IQZ(3)
      MM=IPAM(26)
      DO 20 K=1,3
      IQM=IQZ(K)
      JU=IQM+1
      IF(JU.LT.MM) JU=MM
      IQ=JU
      IQZ(K)=IQ-1
      MMZ(K)=MM
      IF(MM.GT.IQ-3) MMZ(K)=IQ-3
      IGZ(K)=2*MMZ(K)+1
  20  CONTINUE
*
C      write(*,*)"IH(6),TMU==",IH(6),TMU
      IF(IH(6).GE.TMU) THEN
*     ---------------
      TM1=IH(6)+IH(7)
      TI1=IH(6)
      TAP1=IH(7)
      KTAP1=0
      ELSE
*     ---------------
      TM1=TMU+KTAP
      TI1=TMU
      KTAP1=KTAP
      KTAP11=1
      ENDIF
*     ---------------
C      write(*,*)"TI1,ISI,TM1==",TI1,ISI,TM1
      TI1=TI1/ISI
      IREF=TI1
      TM1=TM1/ISI
      IE1=LENGTH-1
      ILT=TM1+IPAM(41)
*     IF(ILT.GT.NE) GO TO 900
C       write(*,*)"ILT,IE1==",ILT,IE1
       IF( ILT .GT.IE1) THEN
       IH(11)=0
       DO 1347 KI=1,IE1
       E(KI)=0.0
 1347  CONTINUE
       RETURN
       ENDIF
      C=0
      IF (IQQ.NE.0) THEN
*     ------------------
      C=IQQ-IH(16)/ISI
      ENDIF
*     ----
*
*++      Defin limits for first gate
*
      I=1
      J=1
      B=TT1(1)+C
      IF(TM1.LE.B) THEN
*     -----------------
      T11(J)=TT1(I)+C
      T21(J)=TT2(I)+C
      ELSE
*     -----------------
      T11(J)=TM1
      T21(J)=TM1+TT2(I)-TT1(I)
      ENDIF
*     -----------------
      INOM=IPAM(41)
      IF(IE1.LT.T21(J)) T21(J)=IE1
      IDIF=T21(J)-T11(J)
      KP=1.5*INOM
      INDEX1(J)=I
      IF(IDIF.LT.KP) THEN
      T21(J)=T11(J)+KP
      ENDIF
*     IF(IDIF.LT.KP) GO TO 900
*
*++       Procesing next gate if any
*
  51  CONTINUE
      J=J+1
  52  CONTINUE
      I=I+1
      IF(I.LE.N) GO TO 60
      NU=J-1
      GO TO 54
  60  CONTINUE
      IDIF=TT1(I)-T11(J-1)+C
      IF(IDIF.LT.INOM) GO TO 52
      T11(J)=TT1(I)+C
      T21(J)=TT2(I)+C
      I1=(T11(J-1)+T21(J-1)+INOM)/2
      I2=(T11(J)+T21(J))/2
      IF(I1.GT.I2) GO TO 52
      INDEX1(J)=I
      IDIF=IE1-T21(J)
      IF(IDIF.LT.0)T21(J)=IE1
      IDIF=T21(J)-T11(J)
      KP=1.5*INOM
      IF(IDIF.GE.KP) GO TO 51
      NU=J-1
  54  CONTINUE
*
*
*++        Computeing start times for application windows
*
*
*    update  time  91.12.17
*
         IBT=IPAM(28)
         IF( IBT .LE. TI1 ) THEN
*        -----------------------
         TI1= TI1-IBT +1
         ELSE
*        ----
         TI1=1
         ENDIF
*        -----
      TT3(1)=TI1
      DO 56 J=1,NU
      IF(J.NE.1) TT3(J)=(T11(J-1)+T21(J-1))/2
      TT4(J)=(T11(J+1)+T21(J+1))/2
      IF(J.EQ.NU)GO TO 59
      GO TO 56
  59  TT4(J)=IE1
  56  CONTINUE
      N=NU
*
*
      IV(3)=NU
      IV(4)=0
      IV(5)=0
      IV(6)=0
      IV(7)=0
      IV(8)=NLE
      IV(9)=0
      IV(10)=0
*
*
      K=1
*    TT3(1)=0
*       TT3(1)=0
      IREF=0
      DO 301 J=1,NU
      K=INDEX1(J)
      LOP(J)=IGZ(K)
      IQ=IQZ(K)+1
      IDEL(J)=IQ-MMZ(K)-1
      IV(100+J)=IDEL(J)
      IV(110+J)=LOP(J)
      IV(90+J)=LOP(J)-1
      IV(10+J)=T11(J)
      IV(20+J)=T21(J)-IV(10+J)+1
      IV(30+J)=INOM
      IV(80+J)=TT4(J)-TT3(J)+1
      IPAM(90+J)=IQZ(K)+1
         IF(J.EQ.1) THEN
*        ------------------
         IV(41)=TT3(1)
         IV(51)=TT4(1)-TT3(1)+1
         ELSE
*        ------------------
         IV(40+J)=TT3(J)-(LOP(J)+IDEL(J))+1
         IV(50+J)=TT4(J)-IV(40+J)+1
         ENDIF
*        ------------------
  301 CONTINUE
      IV(91)=-IDEL(1)+1
*     IV(61)=0
       IV(61)=TT3(1)
      IV(71)=0
      IF(NU.NE.1) THEN
*     ----------------
      DO 302 J=2,NU
      IV(60+J)=TT3(J)
      IV(70+J)=TT4(J-1)-TT3(J)+1
  302 CONTINUE
      ENDIF
*     ---------------
*
*
      DO 303 J=1,NU
       K=INDEX1(J)
      IV(120+J)=IPAM(15+K)
  303 CONTINUE
*
*
*
      J=0
*
      IV(3)=NU
*
C       IF(IPAM(32).EQ.1) THEN
*      ----------------------
C        write(*,*)"Computed parameters :",(IV(IJ),IJ=1,180)
C       ENDIF
*      ----------------------
*
        DO 457 K=1,NE1
        TRF(K)=E(K)
  457   CONTINUE
*
C       write(*,*)"TRF= :",(TRF(IJ),IJ=1001,1010)
*
* @brief 完成自相关、解Toeplitz方程和反褶积等功能。
* @param[in] IV  工作缓冲区，存放计算参数
* @param[in] TRF  工作缓冲区，存放输入地震道
* @param[out] ZS  工作缓冲区，存放反褶积结果
* @param[] ZSC  工作缓冲区，存放中间结果
**
       CALL PRDECON_CORE_DECON (IV,TRF,ZS,ZSC)
c       write(*,*)"ZS= :",(ZS(IJ),IJ=1001,1010)
c       write(*,*)"NLE,NE1=",NLE,NE1
*
*
      IF(NLE.LT.NE1) THEN
         DO 2 I=NLE,NE1
            ZS(I)=0.
    2    CONTINUE
      ENDIF
*
      RETURN
*
      END

*
*
* @brief 通过插值计算算子参数。
* @param[in] INHED  当前空变点值
* @param[in] IFUNC  空变组数
* @param[in] IWIN   时窗个数
* @param[in/out] IPAM   公共参数表
* @param[in] IPAR   空变参数
*
      SUBROUTINE PRDECON_CORE_CALPARM(INHED,IFUNC,IWIN,IPAM,IPAR)
      DIMENSION IPAM(*),IPAR(*)
      N=(IFUNC-1)*25
      K=0
      IF(INHED.LE.IPAR(1)) THEN
      K=1
      GOTO 10
      ENDIF
      IF(INHED.GE.IPAR(N+2)) THEN
      K=IFUNC
      GOTO 10
      ENDIF
      DO I=1,IFUNC-1
      N1=(I-1)*25
      N2=I*25
      IF(INHED.GE.IPAR(N1+1).AND.INHED.LE.IPAR(N1+2)) THEN
       K=I
       GOTO 10
      ENDIF
      IF(INHED.GE.IPAR(N1+2).AND.INHED.LT.IPAR(N2+1)) THEN 
       K1=I
       K2=I+1
       GOTO 20
      ENDIF
      ENDDO
      K=IFUNC
10    CONTINUE
      N=(k-1)*25+5
      DO I=1,20
      IPAM(i)=IPAR(N+I)
      ENDDO
      RETURN
20    CONTINUE
      n1=(k1-1)*25
      n2=(k2-1)*25
      iph1=ipar(n1+2)
      iph2=ipar(n2+1)
      DO I=1,IWIN
      n1=n1+1
      n2=n2+1
      iw1=ipar(n1+5)+(ipar(n2+5)-ipar(n1+5))/((iph2-iph1)*1.0)*
     1     (INHED-iph1)
      iw2=ipar(n1+10)+(ipar(n2+10)-ipar(n1+10))/((iph2-iph1)*1.0)*
     1     (INHED-iph1)
      PD=ipar(n1+15)+(ipar(n2+15)-ipar(n1+15))/((iph2-iph1)*1.0)*
     1     (INHED-iph1)
      cf=ipar(n1+20)+(ipar(n2+20)-ipar(n1+20))/((iph2-iph1)*1.0)*
     1     (INHED-iph1)
      IPAM(I)=IW1
      IPAM(I+5)=IW2
      IPAM(I+10)=PD
      IPAM(I+15)=CF
      ENDDO
      RETURN
      END

*
* @brief 完成自相关、解Toeplitz方程和反褶积等功能。
* @param[in] IV  工作缓冲区，存放计算参数
* @param[in] TRF  工作缓冲区，存放输入地震道
* @param[out] ZS  工作缓冲区，存放反褶积结果
* @param[] ZSC  工作缓冲区，存放中间结果
*
      SUBROUTINE PRDECON_CORE_DECON (IV,TRF,ZS,ZSC)
      DIMENSION  IV(*),TRF(*),ZS(*),ZSC(*)
      DIMENSION  COEF(8000),F(10),B(500)

      DO 10 I=1,5
         F(I)=IV(120+I)/1000.
   10 CONTINUE
*
          IST=IV(61)+10
         DO 12 I=1,IST
         ZS(I)=0.0
 12      CONTINUE
*
*
* @brief     NU是算子个数
*
      NU=IV(3)

      DO  20  I = 1 , NU
        L=IV(30+I)
        IW1=IV(10+I)+1
        IW2=IV(20+I)
C        write(*,*)"L,IW1,IW2=",L,IW1,IW2
*
* @brief 完成自相关功能。
* @param[out] COEF  工作缓冲区，存放自相关结果
* @param[in] L  COEF缓冲区长度
* @param[in] TRF（IW1）  工作缓冲区，存放输入地震道
* @param[in] IW1  时窗起始地址
* @param[in] IW2  时窗长度
*
c      write(*,*) "IW1,IW2,L==",IW1,IW2,L
      CALL PRDECON_CORE_DECON_COR(COEF,L,TRF(IW1),IW2,TRF(IW1),IW2,1,0)
*
c         write(*,*) "Autocor function:",(COEF(IJ),IJ=1,L)
         IF(COEF(1).LT.0.001) THEN
         COEF(1)=1.0
         ENDIF
*
        CX=COEF(1)
        DO 15 J=1,L
           COEF(J)=COEF(J)/CX
   15   CONTINUE
*
c        write(*,*) "F(I)==", F(I)
        COEF(1)=COEF(1)*F(I)
*
        L1=IV(110+I)
        IGG=IV(100+I)
c        write(*,*) "L1,IGG=",L1,IGG
        DO 16 J=1,L1
        B(J)=COEF(IGG+J)
   16   CONTINUE
*
* @brief 完成解Toeplitz方程功能。
* @param[out] ZSC  工作缓冲区，存放算子
* @param[in] COEF  工作缓冲区，存放自相关结果
* @param[in] B  工作缓冲区，存放延迟预测步长的自相关结果
* @param[in] L1  算子长度
* @param[out] IRE  返回标志
*
c         write(*,*) "L1==",L1
         CALL PRDECON_CORE_DECON_TOBLIZE (ZSC,COEF,B,L1,IRE)
*
C        IF( IPAM(32).EQ.1)THEN
c          write(*,*) "Predictive operator:",(zsc(ij),ij=1,l1)
C        ENDIF
*
           DO 3009 K=1,L1
           ZSC(K)=-1*ZSC(K)
 3009      CONTINUE

            DO 902 K=1,L1
            COEF(K)=0.
 902        CONTINUE
        L1=IV(80+I)
        L2=IV(50+I)
        L3=IV(110+I)
        L4=IV(90+I)
        IW1=IV(40+I)+1
         IF(L4.LT.0) THEN
*        -----------------
         L4=0
         ENDIF
*        ------
*
c       write(*,*)"IW1,L1,L2,L3,L4=",IW1,L1,L2,L3,L4
*
*
* @brief 完成褶积功能。
* @param[out] COEF  工作缓冲区，存放褶积结果
* @param[in] L1  COEF缓冲区长度
* @param[in] TRF（IW1）  工作缓冲区，存放输入地震道
* @param[in] IW1  时窗起始地址
* @param[in] L2  时窗长度
* @param[in] ZSC  工作缓冲区，存放预测算子
* @param[in] L3  预测算子长度
* @param[in] L4  算子位移量
*
      CALL PRDECON_CORE_DECON_COR ( COEF,L1,TRF(IW1),L2,ZSC,L3,-1,L4 )
C           IF( IPAM(32).EQ.1) THEN
C           KO1=L1-100
C           ENDIF
*
          IF(IV(90+I).LT.0) THEN
*         ----------------------
          L4=-IV(90+I)
          LL=L1-L4-1
          DO 903 K=1,LL
          COEF(L1-K+1)=COEF(LL-K+1)
 903       CONTINUE
*
          L5=-IV(90+I)+1
           DO 904 K=1,L5
           COEF(K)=0.
 904       CONTINUE
*
           ENDIF
*          --------
*
C        write(*,*) "I==",I
        IF (I.NE.1) THEN
*       ----------------
           LINT=IV(70+I)
           INT1=IV(60+I)+1
*
*
* @brief 过渡区斜坡加权处理。
* @param[in]  ZS      工作缓冲区一
* @param[in/out]  COEF  工作缓冲区二
* @param[in]  LINT      过渡区长度
* @param[in]  INT1   X缓冲区起始位置
*
           CALL PRDECON_CORE_DECON_INT ( ZS,COEF,LINT,INT1 )
        ENDIF
*       -----
        L1=IV(80+I)
        IW1=IV(60+I)+1
*
C        write(*,*) "IW1,L1==",IW1,L1
        DO 19 J=1,L1
           ZS(IW1+J-1)=COEF(J)
   19   CONTINUE
*
   20 CONTINUE
*
      RETURN
      END
*
* @brief 完成相关或褶积功能。
* @param[out] Y  工作缓冲区，存放相关或褶积结果
* @param[in] LY  Y缓冲区长度
* @param[in]  X  工作缓冲区，存放输入地震道
* @param[in] LX  X缓冲区长度
* @param[in] F   工作缓冲区，存放预测算子
* @param[in] LF  X缓冲区长度
* @param[in] IFLAG  相关或褶积标志。 >=0 ：相关；其它：褶积
* @param[in] ISHIFT  算子位移量
*
      SUBROUTINE PRDECON_CORE_DECON_COR (Y,LY,X,LX,F,LF,IFLAG,ISHIFT)
      DIMENSION  Y(*),X(*),F(*)

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

      DO  100  I = 1 , LY
*
          SUM1=0.0
          IS=1
          IFL=LF
          IA=II

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
          IF (IS.GT.IFL) STOP 9999
*
          DO 50 J=IS,IFL
             IA=IA+1
             SUM1=SUM1+X(IA)*F(J)
   50     CONTINUE
          II=II+1
          Y(I)=SUM1
  100 CONTINUE
      RETURN
      END
*
* @brief 过渡区斜坡加权处理。
* @param[in]  X      工作缓冲区一
* @param[in/out]  Y  工作缓冲区二
* @param[in]  L      过渡区长度
* @param[in]  INT1   X缓冲区起始位置
*
*:
*:         _           _
*:         |\         /|
*:         | \       / |
*:         |  \     /  |
*:         |   \   /   |
*:         |    \ /    |
*:         |     X     |
*:         |    / \    |
*:         |   /   \   |
*:         |  /     \  |
*:         | /       \ |
*:         |/         \|
*:         +-----------+
*:
      SUBROUTINE PRDECON_CORE_DECON_INT (X,Y,L,INT1)
      DIMENSION  X(*),Y(*)
      IF(L.LE.1) RETURN
      XL=L-1.0
      XL=1.0/XL
      DO 10 I=0,L-1
         Y(I+1)=X(INT1+I)+(Y(I+1)-X(INT1+I))*I*XL
   10 CONTINUE
      RETURN
      END
*
*
*
*
*===============================================================*
*:                                                              *
*: 功能：                                                        *
*: --------                                                     *
*:    解Toeplitz方程                                             *
*:                                                              *
*===============================================================*
*
* @brief 完成解Toeplitz方程功能。
* @param[out] X  工作缓冲区，存放算子
* @param[in] R  工作缓冲区，存放自相关结果
* @param[in] G  工作缓冲区，存放延迟预测步长的自相关结果
* @param[in] LC  算子长度
* @param[out] IER  返回标志
*
*
      SUBROUTINE PRDECON_CORE_DECON_TOBLIZE(X,R,G,LC,IER)
      DIMENSION X(*),R(*),G(*)
      DIMENSION A(2048),B(2048)
*
      IER=1
      DO 10 I=1,LC
      A(I)=0.
      B(I)=0.
      X(I)=0.
 10     CONTINUE
      IF (ABS(R(1)).LT.0.000001) THEN
*      -------------------------------
      IER=0
      RETURN
      ENDIF
*      -----
      A(1)=1.0/R(1)
      X(1)=G(1)*A(1)
      DO 1 I=1,LC-1
      R0=0.
      T0=0.
      DO 20 K=1,LC
      B(K)=A(K)
 20     CONTINUE
      DO 2 J=1,I
      R0=R0+R(2+I-J)*B(J)
      T0=T0+R(2+I-J)*X(J)
 2      CONTINUE
      Q=1.0/(1.0-R0*R0)
      DO 3 J=1,I+1
      A(J)=Q*(B(J)-R0*B(I+2-J))
 3      CONTINUE
      DO 4 J=1,I+1
      X(J)=X(J)+(G(I+1)-T0)*A(I+2-J)
 4      CONTINUE
 1      CONTINUE
      RETURN
      END
*
*
* @brief 根据选件是否做振幅调整。
* @param[in] IPAM 公共参数表
* @param[in/out] E  地震道
* @param[] ZS    工作缓冲区
* @param[] IZON  工作缓冲区
* @param[] IZOM  工作缓冲区
* @param[in/out] IH  道头
*
*
      SUBROUTINE PRDECON_ADJUST(IPAM,E,ZS,IZON,IZOM,IH)
      DIMENSION IPAM(*),E(*),ZS(*),IH(*)
      REAL IZON(*),IZOM(*)
      EQUIVALENCE (FM,IFM)
*
*
      IF( IH(11).NE.1) RETURN
      ISI=IH(9)
      IF(IPAM(37).EQ.1)THEN
*--------------------------
      KTZ=IPAM(23)
         IF(IH(19).NE.0) THEN
*        ----------------------
         KTZ=KTZ+IPAM(161)
         ENDIF
*        ----------------------
      DO 10 I=1,KTZ
      ZS(I)=0
  10  CONTINUE
      NTAP=200./ISI
      FZ=1./NTAP
      F=0.
      KI=KTZ+1
      KF=KTZ+NTAP
      DO 20 I=KI,KF
      F=F+FZ
      ZS(I)=ZS(I)*F
  20  CONTINUE
      ENDIF
*     --------------
      NE1=IH(10)/IH(9)
*
* @brief    IF(IPAM(27).EQ.0) 调整振幅
* @brief       否则，直接求和
*
C      write(*,*)"input trace :e=",(e(ij),ij=1001,1010)
C      write(*,*)"filter trace :zs=",(zs(ij),ij=1001,1010)
C      write(*,*)"ipam(27)=",ipam(27)
      IF(IPAM(27).EQ.0) THEN
*     ----------------------
      GX=3.
      NO=IPAM(24)
      INU=IPAM(25)
*
* @brief 振幅调整。
* @param[in] E      地震道
* @param[in] ZS     预测误差道
* @param[in] NE1    道长
* @param[in] NO     调节振幅的时窗长度
* @param[in] INU    时窗位移点数
* @param[in] GX     3.
* @param[out] IZON  输出处理后地震道
* @param[]  IZOM    工作缓冲区
*
      CALL PRDECON_ADJUST_PRO(E,ZS,1,NE1,NO,INU,GX,IZON,IZOM)
*
      DO 30 I=1,NE1
      E(I)=IZON(I)
  30  CONTINUE
      ELSE
*     ----------------------
      DO 40 I=1,NE1
      E(I)=ZS(I)+E(I)
  40  CONTINUE
      ENDIF
*     ----------------------

c       write(*,*)"NE1==",NE1
c       write(*,*)"Result trace :e=",(e(ij),ij=700,750)

*
* @brief 计算数组绝对振幅最大值。
* @param[in]  E    数组
* @param[out] FM   绝对振幅最大值
* @param[in]  NE1  数组长度
*
      CALL PRDECON_ADJUST_MAX(E,FM,NE1)
      IH(14)=IFM
      RETURN
      END

*
*
* @brief 振幅调整。
* @param[in] B2    地震道
* @param[in] B5    预测误差道
* @param[in] IFC   1
* @param[in] JFC   道长
* @param[in] NO    调节振幅的时窗长度
* @param[in] INU   时窗位移点数
* @param[in] GAMAX     3.
* @param[out] XTRITON  输出处理后地震道
* @param[] ZONE    工作缓冲区
*
      SUBROUTINE PRDECON_ADJUST_PRO(B2,B5,IFC,JFC,NO,INU,GAMAX,
     +XTRITON,ZONE)
      IMPLICIT INTEGER*4(C-W)
      REAL GA,GAM,QOC,GAMAX
      DIMENSION B2(*),B5(*),XTRITON(*),ZONE(*)
      DIMENSION XIP(32),XIQ(31)

*      write(*,*) "INU,NO==",INU,NO
      INUNU=2*INU+1
      NOO=2*NO+1
      INUP=INU+1
*
*
* FIRST WINDOW PROCESSING
*
*
      DO 1 I=1,INUNU
        XIP(I)=0
 0001   XIQ(I)=0
*
      XIBB = 0.0
      IA=IFC+INU-1
      NQZ=INUNU*NOO
      NQ=NQZ
*
*      write(*,*) "INUNU,NOO,NQ==",INUNU,NOO,NQ
      DO 2 I=1,NOO
         IA=IA+1
         NP=I-NOO
         IB=IA-INUP
*
         DO 3 J=1,INUNU
            IB=IB+1
*
            XIAB=B2(IA)*B5(IB)
            XIP(J)=XIP(J)+XIAB
            NP=NP+NOO
*            write(*,*) "1J,NP==",J,NP
            ZONE(NP)=XIAB
            XIBB=B5(IB)*B5(IB)
            XIQ(J)=XIQ(J)+XIBB
 0003   CONTINUE
         NQ=NQ+1
*         write(*,*) "1I,NQ==",I,NQ
         ZONE(NQ)=XIBB
 0002 CONTINUE
*
      ICAL=INUP
      XIGAM=0.
*
      DO 5 J=1,INUNU
        IF(XIP(J)-XIGAM) 6,6,5
 0006   XIGAM=XIP(J)
        ICAL=J
 0005 CONTINUE
*
      IF(XIQ(ICAL)) 21,21,22
 0021 GA=0.
      GO TO 9
 0022 GAM=XIGAM
      QOC=XIQ(ICAL)
      GA=-GAM/QOC
*
      IF(GA-GAMAX) 9,9,8
 0008 GA=GAMAX
 0009 IFD=IFC+NO+INU
*
* GENERAL PROCESSING
*
      DO 10 IK=IFC,IFD
*
         XTRITON(IK)=B2(IK)+GA*B5(IK)
*
 0010 CONTINUE
*
      IK=IFD
      IKF=JFC-NO-INU
      IFZ=NOO
 1000 NQ=NQZ
*
*      write(*,*) "IK,IKF,IFZ,NQ==",IK,IKF,IFZ,NQ
      DO 1001 I=1,IFZ
         IK=IK+1
         ICAL=INUP
         XIGAM=0.
         IA=IA+1
         NP=I-NOO
         IB=IA-INUP
*
         DO 1002 J=1,INUNU
            IB=IB+1
*
            XIAB=B2(IA)*B5(IB)
            NP=NP+NOO
*            write(*,*) "2J,NP==",J,NP
            XIP(J)=XIP(J)-ZONE(NP)+XIAB
            ZONE(NP)=XIAB
            IF(XIP(J)-XIGAM) 1006,1006,1007
 1006       XIGAM=XIP(J)
            ICAL=J

1007  CONTINUE
*
         IF(J.EQ.1)GOTO1002
            XIQ(J-1)=XIQ(J)
 1002    CONTINUE
         XIBB=B5(IB)*B5(IB)
         NQ=NQ+1
*         write(*,*) "2I,NQ==",I,NQ
         IF(INUNU.EQ.1)THEN
                           XIQ(1)=XIQ(1)-ZONE(NQ)+XIBB
                        ELSE
                            XIQ(INUNU)=XIQ(INUNU-1)-ZONE(NQ)+XIBB
         ENDIF
         ZONE(NQ)=XIBB
         IF(XIQ(ICAL)) 1021,1021,1022
 1021    GA=0.
         GO TO 1009
 1022    GAM=XIGAM
         QOC=XIQ(ICAL)
         GA=-GAM/QOC
         IF(GA-GAMAX) 1009,1009,1008
 1008    GA=GAMAX
 1009    IL=IK+ICAL-INUP
*
         XTRITON(IK)=B2(IK)+GA*B5(IL)
*
 1001 CONTINUE
*
      IF((IK+NOO)-IKF) 1000,1000,1100
 1100 IF(IK-IKF) 1200,2000,1200
 1200 IFZ=IKF-IK
      GO TO 1000
*
 2000 IK=IK+1
*
* @brief 处理最后一个时窗
*
      DO 2001 I=IK,JFC
      XTRITON(I)=B2(I)+GA*B5(I)
*
 2001 CONTINUE
*
      RETURN
      END
*
* @brief 计算数组绝对振幅最大值。
* @param[in]  X    数组
* @param[out] FX   绝对振幅最大值
* @param[in]  L    数组长度
*
      SUBROUTINE PRDECON_ADJUST_MAX (X,FX,L)
      DIMENSION  X(*)
      FX=0.0
      DO 1 I=1,L
         FA=ABS(X(I))
         IF(FX.LT.FA) FX=FA
    1 CONTINUE
      RETURN
      END
