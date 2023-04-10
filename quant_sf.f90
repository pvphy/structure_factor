module array
    implicit none
    double precision, allocatable::mi(:),th(:),ph(:),rwork(:),evl_s(:)
    integer,dimension(:),allocatable::x(:),y(:),x_cl,y_cl
    complex*16,allocatable:: h(:,:),bb(:,:),work(:)
    double precision::nsum_t(1600),nn1(1600),ee1(1600),e_2(1600)
    double precision, allocatable::mi_cl(:),th_cl(:),ph_cl(:)
    doubleprecision::dos_sum_a(1600),dos_sum_b(1600),dos_sum_c(1600),dist1(1600),dist_mi_sum(1600)
end module array

module global
    implicit none
    double precision::t,ua,ub,uc,mu,m_sum,e_sum,u,esum_t,m1,ph1,th1,nsum,e,s_f,sf_A,sf_B_C,m_a,m_b,m_c,sz_i_sz_j,temp
    double precision::lambda,n_a,n_b,n_c,n_a_up,n_b_up,n_c_up,n_a_dn,n_b_dn,n_c_dn,filling_cl,sx_i_sx_j,sy_i_sy_j,sf
    integer::d,nos,ie,flag_isoc,unit_cells,sw,nms
    integer::unit_cells_cl,nos_cl,dim_cl,d_cl
endmodule global

program liebmain
    use array
    use global
    implicit none
    integer::seed,it,i,j
    double precision::temp1,get_mu_s
    double precision::m_temp1,th_temp1,ph_temp1,inpt,filling,fill_inpt
    doubleprecision::avg_qsf_A,avg_qsf_B_C,avg_sf
    !doubleprecision::total_occ_up,total_occ_dn,total_occ,avg_m_a,avg_m_b,avg_m_c,avg_local_m

    open(8,file='input.dat',status='unknown')

    do i=1,8
        read(8,*)inpt
        if (i.eq.1) seed=int(inpt)
        if (i.eq.2) d=int(inpt)
        if (i.eq.3) t=dble(inpt)
        if (i.eq.4) ua=dble(inpt)
        if (i.eq.5) ub=dble(inpt)
        if (i.eq.6) uc=dble(inpt)
        if (i.eq.7) nms=int(inpt)
        if (i.eq.8) fill_inpt=dble(inpt)       
    end do
    close(8)
    unit_cells=(d)**2
    nos=unit_cells
    filling=fill_inpt*2*nos    
      
    print*,'unit cells_____________=',unit_cells
    print*,'sites system___________=',nos
    print*,'hopping parameter______=',t
    print*,'Ua_____________________=',ua
    print*,'ub_____________________=',ub
    print*,'uc_____________________=',uc
    print*,'Monte-Carlo steps______=',nms
    print*,'filling system_________=',filling

    allocate(x(nos),y(nos),h((2*nos),(2*nos)),mi(nos),th(nos),ph(nos))
    
    call lattice_labeling  
    Temp=0.650d0

    
    do iT=1,30
        open(1000+it)   
        if(iT.le.11) Temp=Temp-0.050d0
        if ((iT.gt.11).and.(iT.lt.21)) Temp=Temp-0.010d0
        if ((iT.gt.21)) Temp=Temp-0.0010d0

        avg_qsf_A=0.0d0
        avg_qsf_B_C=0.0d0
        avg_sf=0.0
        do sw=1,nms/10
            j=0
            do i=(sw*nos-nos+1),sw*nos
                j=j+1
                read(1000+it,*) Temp1,m_temp1,th_temp1,ph_temp1 
                mi(j)=m_temp1
                th(j)=th_temp1
                ph(j)=ph_temp1
            enddo

            call matgen
            call diagonalization(h,2,2*nos)
            mu=get_mu_s(filling,temp,2*nos)

            call quan_struc_factor

            avg_sf=avg_sf+(sf/nos**2)    !pi,pi
            deallocate(evl_s)
            
        enddo 

        write(95,*) temp,avg_sf/(nms/10)
        flush(95)
     
    enddo   

endprogram liebmain




subroutine diagonalization(h_temp,flag_diag,dim_1)
    use array
    use global
    implicit none
    integer::lda,lwmax,info,lwork,flag_diag,dim_1,i
    complex*16::h_temp(dim_1,dim_1)

    allocate(rwork(3*(dim_1)-2),work(2*(dim_1)-1),evl_s(dim_1))
    lda=(dim_1)
    lwmax=(dim_1)
    lwork=(2*(dim_1)-1)
    evl_s=0.0d0

    if(flag_diag==1)then
       call zheev('n','u',dim_1,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(flag_diag==2)then
       call zheev('v','u',dim_1,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(info.ne.0)then
        print*,'algorithm failed'  
    endif 
   
    ! do i=1,dim_1
    !     write(726,*) i, evl_s(i)
    ! enddo
    ! stop
    deallocate(rwork,work)  
endsubroutine diagonalization



double precision function get_mu_s(fill,temp2,dim_1)
    use array
    use global
    implicit none
    double precision::f,fL2,fR,mR,mL,m_d,temp2,fill
    integer::i,dim_1
    !integer::fill

    mR = maxval(evl_s)       !right-side chemical potential
    fr=0.0d0
    do i=1,dim_1
        fr=fr+(1.0d0/(exp((evl_s(i)-mR)/Temp2)+1.0d0))
    end do
    mL = minval(evl_s)       !left-side chemical potential
    fL2=0.0d0
    do i=1,dim_1
        fL2=fL2+(1.0d0/(exp((evl_s(i)-mL)/Temp2)+1.0d0))
    end do
    m_d = 0.5d0*(mL+mR)    !middle chemical potential
    f=0.0d0
    do i=1,dim_1
        f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
    end do
    !print*,f,fill
    do while(abs(f-fill).ge.1e-8)
        m_d = 0.5d0*(mL+mR)
        f=0.0d0
        do i=1,dim_1
            f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
        end do
        if(f.gt.fill)then
            !if middle filling is above target, make it the new right bound.
            mR = m_d
            fR = f
        elseif(f.lt.fill)then
            !if middle filling is below target, make it the new left bound.
            mL = m_d
            fR = f
        endif
    enddo

    !Return the middle value
    get_mu_s = m_d
    return
end function get_mu_s
	   


subroutine lattice_labeling
    use global
    use array
    implicit none
    integer::ix,iy,i1,i
    double precision::Pi
    pi=4.0*atan(1.0)
    i1=0
    do iy=1,d
        do ix=1,d
            i1=i1+1
            x(i1)=ix                     
            y(i1)=iy
        enddo
    enddo
    i1=0
    do i=1,nos           

        write(12,*) i,x(i),y(i)
        flush(12)
  
    enddo

endsubroutine lattice_labeling

subroutine matgen
    use array
    use global
    implicit none
    integer::l,k,xi,yi,xd,yd,a,b

    h=complex(0.0d0,0.0d0)
    do l=1,nos
            xi=1
            yi=1
            xd=1
            yd=1
            if(x(l).eq.d) xi=-d+1
            if(y(l).eq.d) yi=-d+1

            if(x(l).eq.1) xd=-d+1
            if(y(l).eq.1) yd=-d+1


            do k=1,nos    
                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)+yi)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)-yd)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
  
               
            enddo
        
    enddo




    call matgen_so
    call matgen_hf_u

    do l=1,2*nos
        do k=1,2*nos
            !if(h(l,k).ne.0) write(34,*) l,k,h(l,k)
            if(h(l,k).ne.conjg(h(k,l))) print*,l,k,'not hermitian'
        enddo
    enddo 

endsubroutine matgen



subroutine matgen_mcmf_u
    use array
    use global
    use mtmod
    implicit none
    integer::l,a,b

    m_sum=0.0d0

    do l=1,nos
        a=2*l-1
        b=2*l-1
        h(a,b)=(u/2.0d0)*(-mi(l)*cos(th(l)))+(u/2.0d0)
        h(a,b+1)=(-u/2.0d0)*mi(l)*sin(th(l))*complex(cos(ph(l)),-sin(ph(l)))
        h(a+1,b)=conjg(h(a,b+1))
        h(a+1,b+1)=(u/2.0d0)*(-mi(l)*(-cos(th(l))))+(u/2.0d0)

        m_sum=m_sum+((mi(l))**2.0d0)*(u/4.0d0)
    enddo

endsubroutine matgen_mcmf_u


subroutine quan_struc_factor
    use array
    use global
    implicit none
    double precision::pi,m_x,m_y,m_z,mimj,qx,qy,rx,ry
    integer::i,j
    pi=4.0*atan(1.0)

    sf=0.0d0
    do i=1,nos ! nos is no of sites  
        do j=1,nos
            qx=pi
            qy=pi
            
            call s_xy_corr(i,j)
            call sz_corr(i,j)
            m_x=sx_i_sx_j
            m_y=sy_i_sy_j
            m_z=sz_i_sz_j
            mimj=m_x+m_y+m_z
            rx=qx*(x(i)-x(j))
            ry=qy*(y(i)-y(j))
            sf=sf+exp(complex(0.0d0,rx+ry))*(mimj)  

        enddo
    enddo
  
endsubroutine quan_struc_factor

subroutine s_xy_corr(i,j)
    use global
    use array
    implicit none
    integer::l,i,j
    doubleprecision::c_i_dn_cdag_i_up,c_j_dn_cdag_j_up,c_i_dn_cdag_j_up,c_i_up_cdag_i_dn,cdag_i_dn_c_j_up 
    double precision::cdag_i_up_c_j_up,cdag_i_dn_c_j_dn,cdag_i_up_c_j_dn,c_i_up_cdag_j_dn,c_j_up_cdag_j_dn
    doubleprecision::c_i_up_cdag_j_up,fermi_fn,sx1,sx2,sx3,sx4

    ! sx_i_sx_j=0.0d0
    ! sy_i_sy_j=0.0d0
    
    !do i=1,nos
        !do j=1,nos
 

           
        
            
            c_i_dn_cdag_i_up=0.0d0
            c_j_dn_cdag_j_up=0.0d0
            c_i_dn_cdag_j_up=0.0d0
            c_i_up_cdag_i_dn=0.0d0
            c_i_up_cdag_j_up=0.0d0
            c_j_up_cdag_j_dn=0.0d0
            c_i_up_cdag_j_dn=0.0d0
            cdag_i_up_c_j_dn=0.0d0
            cdag_i_dn_c_j_dn=0.0d0
            cdag_i_up_c_j_up=0.0d0
            cdag_i_dn_c_j_up=0.0d0


            do l=1,2*nos

                fermi_fn=(1/(1.0d0+(exp((evl_s(l)-mu)/temp))))

                c_i_dn_cdag_i_up=c_i_dn_cdag_i_up+((h(2*i,l))*conjg(h(2*i-1,l)))*(1-fermi_fn)   !1

                c_j_dn_cdag_j_up=c_j_dn_cdag_j_up+((h(2*j,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !2

                c_i_dn_cdag_j_up=c_i_dn_cdag_j_up+((h(2*i,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !3

                c_i_up_cdag_i_dn=c_i_up_cdag_i_dn+((h(2*i-1,l))*conjg(h(2*i,l)))*(1-fermi_fn)    !4

                c_i_up_cdag_j_up=c_i_up_cdag_j_up+((h(2*i-1,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)   !6

                c_j_up_cdag_j_dn=c_j_up_cdag_j_dn+((h(2*j-1,l))*conjg(h(2*j,l)))*(1-fermi_fn)     !9

                c_i_up_cdag_j_dn=c_i_up_cdag_j_dn+((h(2*i-1,l))*conjg(h(2*j,l)))*(1-fermi_fn)     !10

                cdag_i_up_c_j_dn=cdag_i_up_c_j_dn+((h(2*j,l))*conjg(h(2*i-1,l)))*fermi_fn      !5

                cdag_i_dn_c_j_dn=cdag_i_dn_c_j_dn+((h(2*j,l))*conjg(h(2*i,l)))*fermi_fn       !7

                cdag_i_up_c_j_up=cdag_i_up_c_j_up+((h(2*j-1,l))*conjg(h(2*i-1,l)))*fermi_fn   !8

                cdag_i_dn_c_j_up=cdag_i_dn_c_j_up+((h(2*j-1,l))*conjg(h(2*i,l)))*fermi_fn     !11

            enddo

        
        
            sx1=(c_i_dn_cdag_i_up*c_j_dn_cdag_j_up)+(c_i_dn_cdag_j_up*cdag_i_up_c_j_dn)

            sx2=(c_i_up_cdag_i_dn*c_j_dn_cdag_j_up)+(c_i_up_cdag_j_up*cdag_i_dn_c_j_dn)

            sx3=(c_i_dn_cdag_i_up*c_j_up_cdag_j_dn)+(c_i_dn_cdag_j_up*cdag_i_up_c_j_up)

            sx4=(c_i_up_cdag_i_dn*c_j_up_cdag_j_dn)+(c_i_up_cdag_j_dn*cdag_i_dn_c_j_up)


    


                    
           
            sx_i_sx_j=(sx1+sx2+sx3+sx4)*(0.250d0)          !*((-1)**(abs(i-j)))
            sy_i_sy_j=(-sx1+sx2+sx3-sx4)*(0.250d0)          !*((-1)**(abs(i-j)))


        !enddo
    !enddo

endsubroutine s_xy_corr

subroutine sz_corr(i,j)
    use array
    use global
    implicit none
    integer::l,i,j
    complex*16::c_i_up_cdag_i_up,c_j_up_cdag_j_up,c_i_dn_cdag_j_up,c_j_dn_cdag_j_dn,cdag_i_dn_c_j_up 
    complex*16::cdag_i_up_c_j_up,c_i_dn_cdag_j_dn,cdag_i_up_c_j_dn,c_i_up_cdag_j_dn,c_i_dn_cdag_i_dn
    complex*16::c_i_up_cdag_j_up,cdag_i_dn_c_j_dn,sz1,sz2,sz3,sz4
    double precision::fermi_fn

   ! sz_i_sz_j=0.0d0
            
    !do i=1,nos
        !do j=1,nos

                    
    
            c_i_up_cdag_i_up=0.0d0
            c_j_up_cdag_j_up=0.0d0
            c_i_up_cdag_j_up=0.0d0
            c_i_dn_cdag_i_dn=0.0d0
            c_j_dn_cdag_j_dn=0.0d0
            c_i_dn_cdag_j_up=0.0d0
            c_i_up_cdag_j_dn=0.0d0
            c_i_dn_cdag_j_dn=0.0d0
            cdag_i_up_c_j_up=0.0d0
            cdag_i_dn_c_j_up=0.0d0
            cdag_i_up_c_j_dn=0.0d0
            cdag_i_dn_c_j_dn=0.0d0   

            do l=1,2*nos
                fermi_fn=(1.0d0/(1.0d0+(exp((evl_s(l)-mu)/temp))))
                !if((i.eq.89).and.(j.eq.89)) write(*,*) fermi_fn,temp,evl_s(l),mu,h(2*i-1,l)

                c_i_up_cdag_i_up=c_i_up_cdag_i_up+((h(2*i-1,l))*conjg(h(2*i-1,l)))*(1-fermi_fn)   !1

                c_j_up_cdag_j_up=c_j_up_cdag_j_up+((h(2*j-1,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !2

                c_i_up_cdag_j_up=c_i_up_cdag_j_up+((h(2*i-1,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)    !3

                c_i_dn_cdag_i_dn=c_i_dn_cdag_i_dn+((h(2*i,l))*conjg(h(2*i,l)))*(1-fermi_fn)    !4

                c_j_dn_cdag_j_dn=c_j_dn_cdag_j_dn+((h(2*j,l))*conjg(h(2*j,l)))*(1-fermi_fn)   !5

                c_i_dn_cdag_j_up=c_i_dn_cdag_j_up+((h(2*i,l))*conjg(h(2*j-1,l)))*(1-fermi_fn)     !6

                c_i_up_cdag_j_dn=c_i_up_cdag_j_dn+((h(2*i-1,l))*conjg(h(2*j,l)))*(1-fermi_fn)     !7

                c_i_dn_cdag_j_dn=c_i_dn_cdag_j_dn+((h(2*i,l))*conjg(h(2*j,l)))*(1-fermi_fn)      !8

                cdag_i_up_c_j_up=cdag_i_up_c_j_up+((h(2*j-1,l))*conjg(h(2*i-1,l)))*fermi_fn       !9

                cdag_i_dn_c_j_up=cdag_i_dn_c_j_up+((h(2*j-1,l))*conjg(h(2*i,l)))*fermi_fn   !10

                cdag_i_up_c_j_dn=cdag_i_up_c_j_dn+((h(2*j,l))*conjg(h(2*i-1,l)))*fermi_fn     !11

                cdag_i_dn_c_j_dn=cdag_i_dn_c_j_dn+((h(2*j,l))*conjg(h(2*i,l)))*fermi_fn   !12
            enddo
            
            sz1=(c_i_up_cdag_i_up*c_j_up_cdag_j_up)+(c_i_up_cdag_j_up*cdag_i_up_c_j_up)
            
            sz2=(c_i_dn_cdag_i_dn*c_j_up_cdag_j_up)+(c_i_dn_cdag_j_up*cdag_i_dn_c_j_up)
            
            sz3=(c_i_up_cdag_i_up*c_j_dn_cdag_j_dn)+(c_i_up_cdag_j_dn*cdag_i_up_c_j_dn)
            
            sz4=(c_i_dn_cdag_i_dn*c_j_dn_cdag_j_dn)+(c_i_dn_cdag_j_dn*cdag_i_dn_c_j_dn)

                        
            
            sz_i_sz_j=real(sz1-sz2-sz3+sz4)*(0.250d0)!*((-1)**(abs(i-j)))
            
        
        !enddo
    !enddo
endsubroutine sz_corr
