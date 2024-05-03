!
module condiciones
                implicit none 
                integer :: N,nt
                double precision :: p0, pf, upot, dt
                double precision, allocatable, dimension(:) :: ut
end module condiciones

module posiciones
        implicit none
        double precision, allocatable, dimension(:) :: rx,ry,rz
end module posiciones

module velocidades
        implicit none
        double precision, allocatable, dimension(:) :: vx,vy,vz
end module velocidades

module fuerzas
        implicit none
        double precision, allocatable, dimension(:) :: fx,fy,fz
end module fuerzas

program practica  
        use condiciones
        use posiciones
        use velocidades
        use fuerzas
        call leerConfiguracion
        allocate(rx(N),ry(N),rz(N))
        allocate(vx(N),vy(N),vz(N))
        allocate(fx(N),fy(N),fz(N))
        call generarPosiciones
        call guardarPosiciones
        call generarVelocidades
        call calcularFuerzas
        call ejecutarSimulacion
        write(*,*) "La suma de fuerzas es: ", sum(fx+fy+fz)
        write(*,*) "La energia potencial es: ", upot/N
end program practica

subroutine leerConfiguracion
        use condiciones
        open(1,file="config.txt",status="old",action="read")
        read(1,*) N
        read(1,*) p0
        read(1,*) pf
        read(1,*) nt
        read(1,*) dt
        close(1)
end subroutine

subroutine generarPosiciones
        use condiciones
        use posiciones
        implicit none
        integer :: i,j,k,m
        m = 0
        do i =  0,int(pf)-1
            do j = 0,int(pf)-1
                do k = 0,int(pf)-1
                     m = m+1
                     rx(m) = dble(i)
                     ry(m) = dble(j)
                     rz(m) = dble(k)
                end do
            end do
        end do 
        
end subroutine

subroutine guardarPosiciones
        use condiciones
        use posiciones
        integer :: i
        i = 0
        open(1,file="pos.xyz",status="replace",action="write")
        write(1,*) N
        write(1,*) 
        do i = 1,N
             write(1,100) "C",rx(i),ry(i),rz(i)  
        end do
       ! close(1)
        100 format(a,3f10.5)
end subroutine 

subroutine generarVelocidades
        use velocidades
        use condiciones
        implicit none
        double precision :: r
        integer :: i
        vx = 0.0d0
        vy = 0.0d0
        vz = 0.0d0
        do i = 1,N
                call random_number(r)
                r = 2.0d0*r-1.0d0
                vx(i) = r
                call random_number(r)
                r = 2.0d0*r-1.0d0
                vy(i) = r
                call random_number(r)
                r = 2.0d0*r-1.0d0
                vz(i) = r
        end do
end subroutine

subroutine calcularFuerzas
        use fuerzas
        use velocidades
        use posiciones
        use condiciones
        implicit none
        integer :: i,j
        double precision ::r,eps,sig,pot,dpot,dx,dy,dz
        sig = 1.0d0
        eps = 1.0d0
        fx = 0.0d0
        fy = 0.0d0
        fz = 0.0d0
        pot = 0.0d0
        r = 0.0d0
        upot = 0.0d0
        dpot = 0.0d0
        dx = 0.0d0
        dy = 0.0d0
        dz = 0.0d0
        do i = 1,n-1
                do j = i+1,n
                        dx = rx(i)-rx(j)
                        dy = ry(i)-ry(j)
                        dz = rz(i)-rz(j)
                        ! utilizando condiciones de la minima imagen para evitar colisiones kbronas
                        if (dx>pf*0.50d0) then
                                dx = dx-pf
                        else 
                        if (dx<-pf*0.050d0) dx = dx+pf
                        end if
                        if (dy>pf*0.50d0) then
                                dy = dy-pf
                        else 
                        if (dy<-pf*0.050d0) dy = dy+pf
                        end if
                        if (dz>pf*0.50d0) then
                                dz = dz-pf
                        else 
                        if (dz<-pf*0.050d0) dz = dz+pf
                        end if
                        r = sqrt(dx**2.0d0+dy**2.0d0+dz**2.0d0)
                        call leenardJhones(r,sig,eps,pot,dpot)
                        upot = pot+upot
                        fx(i) = fx(i) + dpot*(dx/r)
                        fy(i) = fy(i) + dpot*(dy/r)
                        fz(i) = fz(i) + dpot*(dz/r)
                        ! calculando las reacciones
                        fx(j) = fx(j)-dpot*(dx/r)
                        fy(j) = fy(j)-dpot*(dy/r)
                        fz(j) = fz(j)-dpot*(dz/r)
                end do
        end do
end subroutine

subroutine leenardJhones(r,sig,eps,pot,dpot)
        implicit none
        double precision, intent(in) :: r,sig,eps
        double precision, intent(out) :: pot,dpot 
        dpot = (24.0d0/r)*(2.0d0*(sig/r)**12.0d0-(sig/r)**6.0d0)
        ! potencial de leenard jhones
        pot = 4.0d0*eps*((sig/r)**12.0d0-(sig/r)**6.0d0)
end subroutine

subroutine ejecutarSimulacion
        use fuerzas
        use velocidades 
        use posiciones
        use condiciones
        integer :: paso,i,f,c
        c = 1
        f = 100
        do paso = 1,nt
                do i = 1,N
                
                        rx(i) = rx(i)+vx(i)*dt + 0.50d0*fx(i)*dt**2.0d0
                        ry(i) = ry(i)+vy(i)*dt + 0.50d0*fy(i)*dt**2.0d0
                        rz(i) = rz(i)+vz(i)*dt + 0.50d0*fz(i)*dt**2.0d0
                        
                        ! condiciones periodicas

                        if(rx(i)>pf)rx(i)=rx(i)-pf
                        if(ry(i)>pf)ry(i)=ry(i)-pf
                        if(rz(i)>pf)rz(i)=rz(i)-pf
                        
                        if(rx(i)<0.0d0)rx(i)=rx(i)+pf
                        if(ry(i)<0.0d0)ry(i)=ry(i)+pf
                        if(rz(i)<0.0d0)rz(i)=rz(i)+pf

                        vx(i) = vx(i) + fx(i)*dt
                        vy(i) = vy(i) + fy(i)*dt
                        vz(i) = vz(i) + fz(i)*dt
                end do
                ukin = 0.50d0*sum(vx**2+vy**2+vz**2)
                ti = 2.0d0*ukin/(3*n)

                
                call calcularFuerzas
               ! if(paso>=f) then
                        write(*,*) paso,upot/dble(N),ti,ukin/dble(N)
                !        c = c+1
                !        f = c*f
                 !       write(*,*) "han pasado ...",f,"..... iteraciones"
                !end if
                write(1,*) N
                write(1,*)
                do i =1,N
                        write(1,100) "C",rx(i),ry(i),rz(i)
                end do
                100 format(a,3f10.5)
        end do
end subroutine
