
      subroutine assign_atoms_to_grid(ndum1,mithrd,mirep,
     &           i1_lic,i2_lic,i3_lic,i4_lic)

      use allocatable_arrays
      implicit real(a-h,o-z)

      integer ndum1 ! how many atoms
      integer mirep ! myrep
      integer i1_lic(ndum1),i2_lic(ndum1),i3_lic(ndum1),i4_lic(ndum1)

      i5=mirep ! 4D
      cll_f_1(:,:,:,:,i5)%num=0

C do all current atom positions first
            
      do iii=1,num_mov_atms

C skip atoms that aren't to be mapped to the grid

        if(i_get_mapped_to_grid(iii).eq.0) cycle

        w1=c_a_x(1,iii,mirep)
        w2=c_a_x(2,iii,mirep)
        w3=c_a_x(3,iii,mirep)
        w4=c_a_x(4,iii,mirep) ! 4D
        if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
          w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
          w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
          w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
          w4=w4-wlen*anint(w4*winv)  ! 4D
        else
          if(w1.lt.xmin.or.w1.gt.xmax.or.
     &       w2.lt.ymin.or.w2.gt.ymax.or.
     &       w3.lt.zmin.or.w3.gt.zmax.or.
     &       w4.lt.wmin.or.w4.gt.wmax) then
            write(*,911)iii,w1,w2,w3,w4
911         format('atom# ',i8,' has coords ',4f10.5,
     &            ' these are out of the box & i_pbc=0; this is fatal')
            stop
          endif
        endif                        ! PBC_NO PBC_YES

C assign the atom to the first grid-pointer

        i1=int((w1-xmin)/x_ff)+1
        i2=int((w2-ymin)/y_ff)+1
        i3=int((w3-zmin)/z_ff)+1
        i4=int((w4-wmin)/w_ff)+1 ! 4D

C make sure that the limits are correctly applied:

        i1=max(i1,num_x_ff)
        i2=max(i2,num_y_ff)
        i3=max(i3,num_z_ff)
        i4=max(i4,num_w_ff)
        i1=min(i1,1)
        i2=min(i2,1)
        i3=min(i3,1)
        i4=min(i4,1)

        i1_lic(iii)=i1
        i2_lic(iii)=i2
        i3_lic(iii)=i3
        i4_lic(iii)=i4 ! 4D
        cll_f_1(i1,i2,i3,i4,i5)%num=cll_f_1(i1,i2,i3,i4,i5)%num+1

      enddo

      do i4=1,num_w_ff
        do i3=1,num_z_ff
          do i2=1,num_y_ff
            do i1=1,num_x_ff
              deallocate(cll_f_1(i1,i2,i3,i4,i5)%ida)
            enddo
          enddo
        enddo
      enddo

      do i4=1,num_w_ff
        do i3=1,num_z_ff
          do i2=1,num_y_ff
            do i1=1,num_x_ff
              allocate(cll_f_1(i1,i2,i3,i4,i5)%ida
     &                (cll_f_1(i1,i2,i3,i4,i5)%num))
            enddo
          enddo
        enddo
      enddo

      cll_f_1(:,:,:,:,i5)%num=0

      do iii=1,num_mov_atms
       cll_f_1(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4_lic(iii),i5)%num=
     & cll_f_1(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4_lic(iii),i5)%num+1
       cll_f_1(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4_lic(iii),i5)%ida(
     & cll_f_1(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4_lic(iii),i5)%num)=
     & iii
      enddo                                           
 
      return
      end

