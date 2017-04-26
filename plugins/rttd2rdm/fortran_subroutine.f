

        SUBROUTINE WRAP(re,im,w,n)

        real(8) :: re(n,n),im(n,n),w(n)
        complex(8) :: c(n,n)
        integer :: p,q,lda,lwork,info
        complex(8) :: work(2*n-1)
        real(8) :: rwork(3*n-2)
        c=0.0d0
        do p=1,n
            do q=1,n
                c(p,q)=cmplx(re(p,q),im(p,q))
            end do
        end do
        !write(*,*)''
        !write(*,*)'Complex:',c
        
        ! call zheev
        lda = n
        lwork = 2*n-1
        info=0

        !write(*,*)'Matrix Re:  ',re
        !write(*,*)'Matrix Im:  ',im
        call zheev('V','U',n,c,lda,w,work,lwork,rwork,info)
        ! unpack eigenvectors
        do p=1,n
            do q=1,n
                re(q,p) = real(c(p,q))
                im(q,p) = real(aimag(c(p,q)))
            end do
        end do

        !write(*,*)'Info:        ',info
        !write(*,*)'Eigenval:    ',w
        !write(*,*)'Eigenvectors:',c
        !write(*,*)'Eigenvec Re: ',re
        !write(*,*)'Eigenvec Im: ',im

        RETURN
        END
