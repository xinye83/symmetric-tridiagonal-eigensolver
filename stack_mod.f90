module stack_mod
  ! a stack data structure for the main program

  implicit none

  type stack
     integer             :: len
     real(8),allocatable :: p(:,:)
     integer,allocatable :: c(:,:)
  end type stack

contains

  subroutine stack_copy(S1,S2)
    ! copy stack S1 to S2

    type(stack) :: S1,S2

    S2%len=S1%len
    allocate(S2%p(S2%len,2),S2%c(S2%len,2))

    S2%p=S1%p
    S2%c=S1%c
  end subroutine stack_copy

  subroutine stack_allo(S)
    ! allocate p and c in stack S
    type(stack) :: S
    allocate(S%p(S%len,2),S%c(S%len,2))
  end subroutine stack_allo

  subroutine stack_deal(S)
    ! deallocate stack S
    type(stack) :: S
    deallocate(S%p,S%c)
  end subroutine stack_deal

  subroutine stack_pop(S,p1,p2,c1,c2)
    ! pop one entry from stack

    type(stack) :: S
    real(8)     :: p1,p2
    integer     :: c1,c2

    type(stack) :: temp

    if (S%len==1) then
       p1=S%p(1,1)
       p2=S%p(1,2)

       c1=S%c(1,1)
       c2=S%c(1,2)

       call stack_deal(S)
       S%len=0
    else  
       call stack_copy(S,temp)
       call stack_deal(S)

       S%len=temp%len-1
       call stack_allo(S)

       S%p=temp%p(2:temp%len,:)
       S%c=temp%c(2:temp%len,:)

       p1=temp%p(1,1)
       p2=temp%p(1,2)

       c1=temp%c(1,1)
       c2=temp%c(1,2)

       call stack_deal(temp)
    end if

  end subroutine stack_pop

  subroutine stack_push(S,p1,p2,c1,c2)
    ! push one entry into stack

    type(stack) :: S
    real(8)     :: p1,p2
    integer     :: c1,c2

    type(stack) :: temp

    if (S%len==0) then
       S%len=1
       call stack_allo(S)
       S%p(1,1)=p1
       S%p(1,2)=p2

       S%c(1,1)=c1
       S%c(1,2)=c2
    else   
       call stack_copy(S,temp)
       call stack_deal(S)

       S%len=temp%len+1
       call stack_allo(S)

       S%p(1,1)=p1
       S%p(1,2)=p2
       S%p(2:S%len,:)=temp%p

       S%c(1,1)=c1
       S%c(1,2)=c2
       S%c(2:S%len,:)=temp%c

       call stack_deal(temp)
    end if

  end subroutine stack_push

  subroutine stack_pick(S,k,p1,p2,c1,c2)
    ! extract k-th entry of stack S
    type(stack) :: S
    integer     :: k,c1,c2
    real(8)     :: p1,p2

    p1=S%p(k,1)
    p2=S%p(k,2)

    c1=S%c(k,1)
    c2=S%c(k,2)

  end subroutine stack_pick

end module stack_mod