program project

  use stack_mod
  use omp_lib
  implicit none
  include 'mpif.h'
  include 'mkl.fi'

  integer,parameter :: master=0
  integer           :: numtasks,rank,ierr,rc,len,stat(MPI_STATUS_SIZE)
  character(MPI_MAX_PROCESSOR_NAME) :: hostname

  character(len=32) :: ifile,ofile
  character(len=10) :: rep
  character(len=7)  :: field
  character(len=19) :: symm

  integer                :: nrows,ncols,nnz
  integer,   allocatable :: indx(:),jndx(:),ival(:)
  real(8),   allocatable :: rval(:)
  complex(8),allocatable :: cval(:)

  real(8)                :: a,b,p1,p2,loc_p,delta_p,temp_p1,temp_p2,res
  real(8)                :: tau,tau2,tau3,eps,tr,lam,tic,toc
  real(8),   allocatable :: alp(:),bet(:),p(:),v(:,:),temp_p(:),loc_e(:),loc_x(:,:)
  real(8),   allocatable :: e(:),x(:,:),temp_x(:,:),loc_temp_e(:),loc_temp_x(:,:)
  real(8),   allocatable :: temp_e(:),orth(:,:)
  integer                :: n,c1,c2,flag,loc_c,delta_c,temp_c1,temp_c2,m,proc
  integer                :: neig,i,j,loc_k,off,maxk
  integer,   allocatable :: c(:),temp_c(:),k(:),order(:)
  type(stack)            :: s1,s2

  character(len=32) :: string
  integer :: ipar(4)
  real(8) :: dpar(2)


  call MPI_INIT(ierr)
  if (ierr .ne. MPI_SUCCESS) then
     write(*,*) 'Error starting MPI program. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  end if

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
  call MPI_GET_PROCESSOR_NAME(hostname, len, ierr)

  !------------------!
  !--- read input ---!
  !------------------!

  call getarg(1,string)
  read(string,*) ipar(1)
  call getarg(2,string)
  read(string,*) ipar(2)
  call getarg(3,string)
  read(string,*) ipar(3)
  call getarg(4,string)
  if (string=='') then
     ipar(4)=0
  else
     ipar(4)=1
     read(string,*) dpar(1)
     call getarg(5,string)
     read(string,*) dpar(2)
  end if

  if (rank==master) then
     write(*,*) '|--------------------------------------|'
     write(*,*) '|--------- input information ----------|'
     write(*,*) '|--------------------------------------|'

     write(*,'(a13,i4)') 'num of nodes',numtasks

     if (ipar(1)==1) then
        write(*,*) 'test matrix: poisson matrix'
     else 
        write(*,*) 'test matrix:'
        write(*,*) '      diagonals are {1,...,100} and sub/super-diagonals are -1'
     end if

     if (ipar(2)/=0) then
        write(*,'(a11,i10)') 'dimension:',ipar(2)
     else
        write(*,*) 'dimension: 10^6 (default)'
     end if

     if (ipar(3)==1) then
        write(*,*) 'eigenvectors will be computed'
     else
        write(*,*) 'eigenvectors will NOT be computed'
     end if

     if (ipar(4)==1) then
        write(*,*) 'search interval is given'
     else
        write(*,*) 'all eigenvalues will be computed'
     end if

     write(*,*) '----------------------------------------'
  end if

  if (ipar(2)/=0) then
     n=ipar(2)
  else
     n=1e6
  end if

  allocate(alp(n),bet(n-1))
  if (ipar(1)==1) then
     alp=2
     bet=-1
  elseif (ipar(1)==2) then
     do i=1,n
        alp(i)=1+(i-1)*99/(n-1)
     end do

     bet=-1
  end if

  ! machine precision
  eps=epsilon(alp(1)/alp(1))

  tau=10*eps   ! tol for cluster in multisection stage
  tau2=5*eps   ! tol for absolute error in eigenvalue
  tau3=eps     ! tol for inverse iteration

  ! ||T||_R
  tr=abs(alp(1))
  do i=2,n
     tr=max(tr,abs(alp(i))+abs(bet(i-1)))
  end do

  if (ipar(4)==1) then
     a=dpar(1)
     b=dpar(2)
  else
     call disk(n,alp,bet,a,b)
  end if

  !---------------------!
  !--- solver starts ---!
  !---------------------!

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  tic=MPI_WTIME()


  if (rank==master) then
     S1%len=0 ! stack for working intervals
     S2%len=0 ! stack for results

     if (ipar(4)==1) then
        call sturm(n,alp,bet,a,c1)
        call sturm(n,alp,bet,b,c2)
     else
        c1=0
        c2=n
     end if

     neig=c2-c1

     write(*,'(a20,e12.3,e12.3)') 'initial interval',a,b
     write(*,'(a20,i12)')       'eigenvalue count',neig

     call stack_push(S1,a,b,c1,c2)

     allocate(p(numtasks),c(numtasks))
  end if

  !--------------------------------!
  !--- isolate each eigenvalues ---!
  !--------------------------------!

  flag=1

  do while (flag)

     if (rank==master) then
        call stack_pop(S1,p1,p2,c1,c2)

        delta_p=(p2-p1)/(numtasks+1)

        do i=1,numtasks
           p(i)=p1+i*delta_p
        end do
     end if

     call MPI_SCATTER(p,1,MPI_REAL8,loc_p,1,MPI_REAL8,master,MPI_COMM_WORLD,ierr)

     call sturm(n,alp,bet,loc_p,loc_c)

     call MPI_GATHER(loc_c,1,MPI_INTEGER,c,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

     if (rank==master) then
        do i=numtasks,0,-1

           if (i==0) then
              temp_p1=p1
              temp_c1=c1
           else
              temp_p1=p(i)
              temp_c1=c(i)
           end if

           if (i/=numtasks) then
              temp_p2=p(i+1)
              temp_c2=c(i+1)
           else
              temp_p2=p2
              temp_c2=c2
           end if

           delta_c=temp_c2-temp_c1

           if ((delta_c>1).AND.(delta_p>tau)) then
              call stack_push(S1,temp_p1,temp_p2,temp_c1,temp_c2)
           else
              if (delta_c>0) then
                 call stack_push(S2,temp_p1,temp_p2,temp_c1,temp_c2)
              end if
           end if

        end do

     end if

     if (rank==master) then
        if (S1%len==0) then
           flag=0

           deallocate(p,c)
        end if
     end if

     call MPI_BCAST(flag,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

  end do

  !---------------------------------------------------------!
  !--- obtain eigenvalue or cluster to desired tolerance ---!
  !---------------------------------------------------------!

  allocate(temp_p(2),temp_c(2))

  if (rank==master) then
     flag=S2%len

     allocate(k(numtasks))
     k=0
     do i=1,flag
        proc=mod(i,numtasks)

        k(proc+1)=k(proc+1)+S2%c(i,2)-S2%c(i,1)
     end do

     maxk=maxval(k)
  end if

  call MPI_BCAST(flag,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(maxk,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_SCATTER(k,1,MPI_INTEGER,loc_k,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

  allocate(loc_e(maxk),loc_x(n,maxk))

  off=0  
  do i=1,flag
     proc=mod(i,numtasks)

     if (rank==master) then
        call stack_pick(S2,i,p1,p2,c1,c2)

        if (proc/=master) then
           temp_p=(/p1,p2/)
           temp_c=(/c1,c2/)

           call MPI_SEND(temp_p,2,MPI_REAL8,proc,1,MPI_COMM_WORLD,ierr)
           call MPI_SEND(temp_c,2,MPI_INTEGER,proc,2,MPI_COMM_WORLD,ierr)
        end if
     end if

     if (rank==proc) then       
        if (rank/=master) then

           call MPI_RECV(temp_p,2,MPI_REAL8,master,1,MPI_COMM_WORLD,stat,ierr)
           call MPI_RECV(temp_c,2,MPI_INTEGER,master,2,MPI_COMM_WORLD,stat,ierr)

           p1=temp_p(1)
           p2=temp_p(2)

           c1=temp_c(1)
           c2=temp_c(2)
        end if

        ! extract eigenvalue(s)
        m=c2-c1
        if (m==1) then
           ! simple eigenvalue
           call bisection(n,alp,bet,p1,p2,c1,c2,tau2,lam)
        else
           lam=(p1+p2)/2
        end if

        ! obtain eigenvector(s) via inverse iteration
        if (ipar(3)==1) then
           allocate(temp_x(n,m))
           call inv_iter(n,alp,bet,lam,m,temp_x,tau3)
        end if

        loc_e(off+1:off+m)=lam

        if (ipar(3)==1) then
           loc_x(:,off+1:off+m)=temp_x
           deallocate(temp_x)
        end if

        off=off+m
     end if

  end do

  if (loc_k/=maxk) then
     loc_e(maxk)=0
     if (ipar(3)==1) then
        loc_x(:,maxk)=0
     end if
  end if

  deallocate(temp_p,temp_c)
  if (rank==master) then
     call stack_deal(S2)
  end if

  !-----------------------------------------!
  !--- gather eigenpairs from each nodes ---!
  !-----------------------------------------!

  if (rank==master) then
     allocate(temp_e(maxk*numtasks),e(neig))
     if (ipar(3)==1) then
        allocate(temp_x(n,maxk*numtasks),x(n,neig))
     end if
  end if

  call MPI_Gather(loc_e,maxk,MPI_REAL8,temp_e,maxk,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
  deallocate(loc_e)

  if (ipar(3)==1) then
     call MPI_Gather(loc_x,maxk*n,MPI_REAL8,temp_x,maxk*n,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
     deallocate(loc_x)
  end if

  if (rank==master) then
     flag=0
     do i=1,numtasks
        e(flag+1:flag+k(i))=temp_e((i-1)*maxk+1:(i-1)*maxk+k(i))

        if (ipar(3)==1) then
           x(:,flag+1:flag+k(i))=temp_x(:,(i-1)*maxk+1:(i-1)*maxk+k(i))
        end if

        flag=flag+k(i)
     end do

     deallocate(temp_e)
     if (ipar(3)==1) then
        deallocate(temp_x)
     end if
  end if

  !----------------------------------------------------------!
  !--- group of close eigenvalues and reorthonomalization ---!
  !----------------------------------------------------------!

!!$  if (ipar(3)==1) then
!!$
!!$     if (rank==master) then
!!$        allocate(order(neig))
!!$        call sort(n,neig,e,order)
!!$     end if
!!$
!!$
!!$  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  toc=MPI_WTIME()

  !-------------------!
  !--- solver ends ---!
  !-------------------!

  ! e: eigenvalues
  ! x: eigenvectors

  !--- computing residuals ---!

  if (rank==master) then
     write(*,'(a13,e10.3)') 'time',toc-tic

     allocate(v(n,neig))
     call symtridmatvec(n,alp,bet,neig,x,v)

     res=0
     do i=1,neig
        v(:,i)=v(:,i)-e(i)*x(:,i)

        res=max(res,dnrm2(n,v(:,i),1))
     end do

     write(*,'(a13,e10.3)') 'max res',res

     allocate(orth(neig,neig))

     call dgemm('T','N',neig,neig,n,DBLE(1),x,n,x,n,DBLE(0),orth,neig)
     do i=1,neig
        orth(i,i)=orth(i,i)-1
     end do

     write(*,'(a13,e10.3)') 'orthog',maxval(abs(orth))
  end if


  call MPI_FINALIZE(ierr)

end program project