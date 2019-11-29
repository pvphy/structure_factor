
!Auxilliary field structure factor
subroutine struc_factor
    use array
    use global
    implicit none
    double precision::pi,m_x,m_y,m_z,mimj,qx,qy,rx,ry
    integer::i,j
    pi=4.0*atan(1.0)
    qx=0.0d0
    qy=0.0d0
    !s_f=0.0d0
    sf_A=0.0d0
    sf_B_C=0.0d0
    do i=1,nos ! nos is no of sites  
        if(tag(i).eq.1)then
            do j=1,nos
                if(tag(j).eq.1)then
                    m_x=(mi(i)*sin(th(i))*cos(ph(i)))*(mi(j)*sin(th(j))*cos(ph(j)))
                    m_y=(mi(i)*sin(th(i))*sin(ph(i)))*(mi(j)*sin(th(j))*sin(ph(j)))
                    m_z=(mi(i)*cos(ph(i)))*(mi(j)*cos(ph(j)))
                    mimj=m_x+m_y+m_z
                    rx=qx*(x(i)-x(j))/2
                    ry=qy*(y(i)-y(j))/2
                    !sf_A=sf_A+(complex(cos(rx+ry),sin(rx+ry)))*(mimj)        
                    sf_A=sf_A+exp(complex(0.0d0,rx+ry))*(mimj)        

                endif
            enddo
        endif

        if((tag(i).eq.2).or.(tag(i).eq.3))then
            do j=1,nos
                if((tag(j).eq.2).or.(tag(j).eq.3))then                  
                    m_x=(mi(i)*sin(th(i))*cos(ph(i)))*(mi(j)*sin(th(j))*cos(ph(j)))
                    m_y=(mi(i)*sin(th(i))*sin(ph(i)))*(mi(j)*sin(th(j))*sin(ph(j)))
                    m_z=(mi(i)*cos(ph(i)))*(mi(j)*cos(ph(j)))
                    mimj=m_x+m_y+m_z
                    rx=qx*(x(i)-x(j))/2
                    ry=qy*(y(i)-y(j))/2
                    !sf_B_C=sf_B_C+(complex(cos(rx+ry),sin(rx+ry)))*(mimj)
                    sf_B_C=sf_B_C+exp(complex(0.0d0,rx+ry))*(mimj)        

                endif
            enddo
        endif
    enddo
endsubroutine struc_factor
