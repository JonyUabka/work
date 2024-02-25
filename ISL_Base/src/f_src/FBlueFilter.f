       module bluefilter_module
        type, public :: bluefilter_obj
        !
         real si ! sample interval(us) 采样间隔
         !integer inited !是否初始化


         integer lf2 ! lf2=((llf/si)/2*2+1)/2
         integer lt  ! length of trace
         integer llf ! filter window length

         REAL, allocatable :: ff(:,:) ! (-lf2:lf2,nw)

         integer nw ! num of windows
         real, allocatable :: time_point(:)
         real, allocatable :: f_low_cut(:)
         real, allocatable :: slope_l(:)
         real, allocatable :: coe(:)
         real, allocatable :: f_high_cut(:)
         real, allocatable :: slope_h(:)
         real, allocatable :: f_expand(:)
         real, allocatable :: f_wide(:)
        end type bluefilter_obj

       !
        type(bluefilter_obj), target :: config(1)

        end module bluefilter_module


        subroutine bluefilter_initialize(nw, llf,
     +    time_point, coe,
     +    f_wide, f_expand,
     +    f_high_cut, slope_h,
     +    f_low_cut, slope_l)
        use bluefilter_module
        implicit none
        integer :: nw, llf
        real :: time_point(nw)
        real :: f_low_cut(nw)
        real :: slope_l(nw)
        real :: coe(nw)
        real :: f_high_cut(nw)
        real :: slope_h(nw)
        real :: f_expand(nw)
        real :: f_wide(nw)
        type(bluefilter_obj), pointer :: p


        p=>config(1)

        p%nw = nw
        p%llf = llf
        !print *,"nw,llf=", nw,llf,time_point
        allocate(p%time_point(p%nw))
        allocate(p%f_low_cut(p%nw))
        allocate(p%slope_l(p%nw))
        allocate(p%coe(p%nw))
        allocate(p%f_high_cut(p%nw))
        allocate(p%slope_h(p%nw))
        allocate(p%f_expand(p%nw))
        allocate(p%f_wide(p%nw))
        p%time_point=time_point
        p%f_low_cut=f_low_cut
        p%slope_l=slope_l
        p%f_high_cut=f_high_cut
        p%slope_h=slope_h
        p%f_expand=f_expand
        p%f_wide=f_wide
        p%coe=coe
        !print *,"nw,llf=", nw,llf,p%time_point
        !print *, "bluefilter_initialize() finished."

        end subroutine bluefilter_initialize


*
*+------------------------------------------------------------
*+    程序功能：对地震道进行蓝色滤波
*+    参数：
*+        trace--     道数据
*+        body_lt--   道头数据长度
*+        nw--        时窗个数
*+        llf--       滤波窗口长度
*+        time_point--时窗值
*+        coe--       加强算子系数
*+        f_wide--    加强频率时窗宽度
*+        f_expand--  加强频率
*+        f_high_cut--高通频率
*+        slope_h--   高通频率斜坡
*+        f_low_cut-- 高通频率
*+        slope_l--   高通频率斜坡
*+------------------------------------------------------------
*+------------------------------------------------------------
      SUBROUTINE bluefilter_process(trace, body_lt,opt)
        use bluefilter_module
        implicit none
        integer body_lt
        real trace(body_lt)
        real opt(10)
        integer it, ifloop, iw
        real taper
        real ifrequency,weight,aa,it1,it2,iw1,iw2

        type(bluefilter_obj), pointer :: p



         integer, pointer :: nw
         real, pointer :: ff(:,:)
         integer, pointer :: lf2, lt
         real, pointer :: si
         real, pointer :: time_point(:),f_low_cut(:)
         real, pointer :: slope_l(:),coe(:)
         real, pointer :: f_high_cut(:),slope_h(:)
         real, pointer :: f_expand(:),f_wide(:)
         REAL :: temp(body_lt),work(body_lt)

        !print *, "bluefilter_process begin"

        p=>config(1)

        ff=>p%ff
        nw=>p%nw
        ! print *,"process nw=",nw,p%nw
        ! time_point=>p%time_point
        lf2=>p%lf2
        lt=>p%lt
        si=>p%si

        lt=opt(5)
        si=opt(6)/1000.

        !print *,"lt,si=", lt,si
        lt=lt/si
        !print *,"lt=", lt
        time_point=>p%time_point
        f_low_cut=>p%f_low_cut
        slope_l=>p%slope_l
        coe=>p%coe
        f_high_cut=>p%f_high_cut
        slope_h=>p%slope_h
        f_expand=>p%f_expand
        f_wide=>p%f_wide

        !if( p%inited .lt. 0) then
        call bluefilter_on_first_data(opt)
        !end if

        !print *, "bluefilter_process opt =" , opt(5)
        !print *, "bluefilter_process lt =", opt(6)
        !print *,"lt=", lt
       do it=1,lt
        temp(it) = trace(it)
        trace(it) = 0.0
       end do
       !print *, "temp(it) = trace(it)"

       it1 = time_point(1)
       !print *, "time_point=",time_point,f_wide

       iw1 = 1

       !print *, "bluefilter_process :: if(nw.eq.1) then"

       if(nw.eq.1) then
        iw2 = lt
       else
        iw2 = time_point(2)
       end if

       do it=1,lt
        work(it) = 0.0
       end do

       do it=iw1,iw2
       do ifloop=-min(lf2,it-1),min(lf2,lt-it)
         work(it) = work(it)+temp(it+ifloop)*p%ff(ifloop,1)
         !print *,"it=",it,iw1,iw2
       end do
       end do
       !print *, "bluefilter_process ::work(it) = work(it)+temp(it+ifloop)*ff(ifloop,1)"

        do it=iw1,it1
        trace(it) = work(it)
        end do

        !print *, "bluefilter_process :: trace(it) = ", work(it)
        if(nw.eq.1) then
        do it=it1+1,iw2
          trace(it) = work(it)
          !print *, "trace1  =" , trace(it)
        end do
        else
        do it=it1+1,iw2
          taper = (iw2-it)/(iw2-it1)
          trace(it) = trace(it)+work(it)*taper
          !print *, "trace2  =" , trace(it)
        end do

        do iw=2,nw

          it1 = time_point(iw)
          iw1 = time_point(iw-1)

          if(iw.eq.nw) then
            iw2 = lt
          else
            iw2 = time_point(iw+1)
          end if

          !print *, "bluefilter_process :: iw2 = time_point(iw+1)"

          do it=1,lt
            work(it) = 0.0
          end do

          do it=iw1,iw2
          do ifloop=-min(lf2,it-1),min(lf2,lt-it)
              work(it) = work(it)+temp(it+ifloop)*p%ff(ifloop,iw)
          end do
          end do

          do it=iw1,it1
            taper =(it-iw1)/(it1-iw1)
            trace(it) = trace(it)+work(it)*taper
            !print *, "trace3  =" , trace(it)
          end do

          if(iw.eq.nw) then
            do it=it1+1,iw2
              trace(it) = work(it)
              !print *, "trace 4  =" , trace(it)
            end do
          else
            do it=it1+1,iw2
              taper = (iw2-it)/(iw2-it1)
              trace(it) = trace(it)+work(it)*taper
              !print *, "trace 5  =" , trace(it)
            end do
          end if

        end do
        end if
        !print *, "bluefilter_process  end "

        do it=1, 50
            !print *, " xxxx end trace  =" , trace(it)
            end do
        end subroutine bluefilter_process

        subroutine bluefilter_finish()
        use bluefilter_module
        IMPLICIT NONE

        type(bluefilter_obj), pointer :: p

         p=>config(1)

         if(allocated(p%time_point)) deallocate(p%time_point)
         if(allocated(p%f_low_cut)) deallocate(p%f_low_cut)
         if(allocated(p%slope_l)) deallocate(p%slope_l)
         if(allocated(p%coe)) deallocate(p%coe)
         if(allocated(p%f_high_cut)) deallocate(p%f_high_cut)
         if(allocated(p%slope_h)) deallocate(p%slope_h)
         if(allocated(p%f_expand)) deallocate(p%f_expand)
         if(allocated(p%f_wide)) deallocate(p%f_wide)
         if(allocated(p%ff)) deallocate(p%ff)

        end subroutine bluefilter_finish



        subroutine bluefilter_on_first_data(opt)
        use bluefilter_module
        implicit none

        real opt(10)
        real lf !todo::onFirstData
        INTEGER llf !todo::onFirstData
        type(bluefilter_obj), pointer :: p

        p=>config(1)
        ! trace length
        !p%lt = opt(5)

        ! sample interval
        !p%si=opt(6)/1000.0

        !print *,"lt si  3=",p%lt,p%si

        !filter window length
        llf= p%llf
        lf = llf/p%si
        lf = lf/2*2+1
        p%lf2   = lf/2

      !
        allocate(p%ff(-p%lf2:p%lf2,p%nw))

       !
       call bluefilter_initMatrix()

       !print *, "bluefilter_on_first_data"

       end subroutine bluefilter_on_first_data


       subroutine bluefilter_initMatrix()
       use bluefilter_module
       implicit none

       integer i,isum,ifloop
       real pi,term,ifn
       integer iw
       real f1,db1,f2,db2,ifc,ifb
       real scl
       real ibpf1,ibpf2,ibpf3,ibpf4,ifc1,ifc2
       real ifrequency,weight,aa,it1,it2,iw1,iw2

       type(bluefilter_obj), pointer :: p

       integer, pointer :: nw
       real, pointer :: ff(:,:)
       integer, pointer :: lf2, lt
       real, pointer :: si
       real, pointer :: time_point(:),f_low_cut(:)
       real, pointer :: slope_l(:),coe(:)
       real, pointer :: f_high_cut(:),slope_h(:)
       real, pointer ::f_expand(:),f_wide(:)

       !print *, "bluefilter_initMatrix"

       p=>config(1)
       ff=>p%ff
       nw=>p%nw
       time_point=>p%time_point
       lf2=>p%lf2
       lt=>p%lt
       si=>p%si

       time_point=>p%time_point
       f_low_cut=>p%f_low_cut
       slope_l=>p%slope_l
       coe=>p%coe
       f_high_cut=>p%f_high_cut
       slope_h=>p%slope_h
       f_expand=>p%f_expand
       f_wide=>p%f_wide

       pi   = 3.1415926
       term = 2*pi*si/1000.
       ifn  = 500./si

       do 1111 iw=1,nw
            f1 = f_low_cut(iw)
            db1 =slope_l(iw)
            f2  = f_high_cut(iw)
            db2 = slope_h(iw)
            ifc = f_expand(iw)
            ifb = f_wide(iw)
            scl =  coe(iw)
            if(f1.gt.f2) then
                print*,'frequency scope overlaping.'
            endif

            if(f1.gt.0.) then
                ibpf1 = f1-f1/db1/2.
            if(ibpf1.le.0) ibpf1=0
                ibpf2 = f1+f1/db1/2.
            if(ibpf2.gt.ifn) ibpf2=ifn
            else
                ibpf1 = 0.
                ibpf2 = 0.
            end if

            if(f2.gt.0.) then
                ibpf3 = f2-f2/db2/2.
            if(ibpf3.lt.ibpf2) ibpf3=ibpf2
            ibpf4 = f2+f2/db2/2.
            if(ibpf4.ge.ifn) ibpf4=ifn
            else
                ibpf3 = ifn-4
                ibpf4 = ifn
            end if

            ifc1 = ifc-ifb
            if(ifc1.lt.ibpf2) ifc1=ibpf2
            ifc2 = ifc+ifb
            if(ifc2.gt.ibpf3) ifc2=ibpf3

            do ifloop=-lf2,lf2
                ff(ifloop,iw) = 0.0
            end do

            if(ibpf1.lt.ibpf2) then
                !isum = ibpf2 - ibpf1 -1
                !do i = 0, isum
                !    ifrequency = ibpf1+i
                do ifrequency=ibpf1,ibpf2-1
                    weight = (ifrequency-ibpf1)/(ibpf2-ibpf1)
                    do ifloop=-lf2,lf2
                        ff(ifloop,iw) = ff(ifloop,iw)+
     +                   cos(term*ifrequency*ifloop)*weight
                    end do
                end do
            end if

        if(ibpf2.lt.ifc1) then
        do ifrequency=ibpf2,ifc1-1
        do ifloop=-lf2,lf2
        ff(ifloop,iw) = ff(ifloop,iw)
     +   +cos(term*ifrequency*ifloop)
        end do
        end do

         end if

        if(ifc1.lt.ifc) then
        do ifrequency=ifc1,ifc-1
            weight = 1+(scl-1)*(ifrequency-ifc1)/(ifc-ifc1)
            do ifloop=-lf2,lf2
                ff(ifloop,iw) = ff(ifloop,iw)+
     +           cos(term*ifrequency*ifloop)*weight
            end do
        end do
        end if

        do ifloop=-lf2,lf2
        ff(ifloop,iw) = ff(ifloop,iw)+cos(term*ifc*ifloop)*scl
        end do

        if(ifc.lt.ifc2) then
        do ifrequency=ifc+1,ifc2
            weight = scl+(1-scl)*(ifrequency-ifc)/(ifc2-ifc)
            do ifloop=-lf2,lf2
                ff(ifloop,iw) = ff(ifloop,iw)+
     +           cos(term*ifrequency*ifloop)*weight
            end do
         end do
         end if

       if(ifc2.lt.ibpf3) then
       do ifrequency=ifc2+1,ibpf3
       do ifloop=-lf2,lf2
        ff(ifloop,iw) = ff(ifloop,iw)
     +   +cos(term*ifrequency*ifloop)
       end do
       end do

       end if

       if(ibpf3.lt.ibpf4) then
        do ifrequency=ibpf3+1,ibpf4
            weight = (ibpf4-ifrequency)/(ibpf4-ibpf3)
            do ifloop=-lf2,lf2
                ff(ifloop,iw) = ff(ifloop,iw)+
     +           cos(term*ifrequency*ifloop)*weight
            end do
       end do

       end if

       aa = ff(0,iw)
       do ifloop=-lf2,lf2
       ff(ifloop,iw) = ff(ifloop,iw)/aa
       end do

1111   continue

       do iw=1,nw
          it1 = time_point(iw)
!          print*,"time_point=",time_point(iw),it1,si
          it1 = it1/si
!          print*,"time_point=",time_point(iw),it1,si
          if(it1.lt.1 ) it1=1
          if(it1.gt.lt) it1=lt
          time_point(iw) = it1

        end do

        if(nw.gt.1) then
        do iw=2,nw
            it1 = time_point(iw-1)
            it2 = time_point(iw)
            !print*,'it1= it2=',it1,it2
            if(it1.ge.it2) then
            print*,'successive controlling point not increase.'
            !pause
            end if
        enddo
        end if
       end subroutine bluefilter_initMatrix


