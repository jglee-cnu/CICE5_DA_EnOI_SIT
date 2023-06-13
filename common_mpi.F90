module common_mpi

      include 'mpif.h'
      integer ierr, mpistatus(mpi_status_size)

contains

      subroutine barrier(label)
       use ice_communicate, only: MPI_COMM_ICE, my_task, master_task
       implicit none
       character(*), optional :: label
       call mpi_barrier(MPI_COMM_ICE, ierr)
       if(my_task==master_task.and.present(label)) print *, '----barrier---',label,'---------'
      end subroutine barrier



      subroutine common_mpi_send(data, target, tag)
       use ice_communicate, only: MPI_COMM_ICE
       implicit none
       real data
       integer target
       integer, optional :: tag
       integer counter, given_tag

       given_tag=0
       if(present(tag)) given_tag=tag
       counter=1

       call mpi_send(data, counter, mpi_real, target, given_tag, MPI_COMM_ICE, ierr)

       if(ierr.ne.0) then
         print *, 'error common_mpi_send count=', counter, 'tag=', given_tag
       endif

      end subroutine common_mpi_send


      subroutine common_mpi_send_1D(data, target, tag)
       use ice_communicate, only: MPI_COMM_ICE
       implicit none
       real data(:)
       integer target
       integer, optional :: tag
       integer counter, given_tag

       given_tag=0
       if(present(tag)) given_tag=tag
       counter=size(data,1)
!       print *, counter, 'counter sended'

       call mpi_send(data, counter*2, mpi_real, target, given_tag, MPI_COMM_ICE, ierr)

       if(ierr.ne.0) then
         print *, 'error common_mpi_send_1D count=', counter, 'tag=', given_tag
         stop
       endif
      
      end subroutine common_mpi_send_1D



      subroutine common_mpi_send_2D(data, target, tag)
       use ice_communicate, only: MPI_COMM_ICE
       implicit none
       real data(:,:)
       integer target
       integer, optional :: tag
       integer counter, given_tag

       given_tag=0
       if(present(tag)) given_tag=tag
       counter=size(data,1)*size(data,2)

       call mpi_send(data, counter, mpi_real, target, given_tag, MPI_COMM_ICE,ierr)

       if(ierr.ne.0) then
           print *, 'error common_mpi_send_2D count=', counter, 'tag=', given_tag
           stop
       endif

      end subroutine common_mpi_send_2D



      subroutine common_mpi_irecv(data, source, tag)
       use ice_communicate, only: MPI_COMM_ICE
       implicit none
       real data
       integer source
       integer, optional :: tag
       integer counter, given_tag
       integer, allocatable, dimension(:) :: request

       given_tag=0
       if(present(tag)) given_tag=tag
       counter=1

       allocate(request(counter))

       call mpi_irecv(data, counter, mpi_real, source, given_tag, MPI_COMM_ICE, request(counter), ierr)

       if(ierr.ne.0) then
         print *, 'error common_mpi_irecv count=', counter, 'tag=', given_tag
       endif

      end subroutine common_mpi_irecv




      subroutine common_mpi_recv(data, source, tag)
       use ice_communicate, only: MPI_COMM_ICE
       implicit none
       real data
       integer source
       integer, optional :: tag
       integer counter, given_tag

       given_tag=0
       if(present(tag)) given_tag=tag
       counter=1

       call mpi_recv(data, counter, mpi_real, source, given_tag, MPI_COMM_ICE, mpistatus, ierr)

       if(ierr.ne.0) then
         print *, 'error common_mpi_recv count=', counter, 'tag=', given_tag
       endif

      end subroutine common_mpi_recv



      subroutine common_mpi_recv_1D(data, source, tag)
       use ice_communicate, only: MPI_COMM_ICE
       implicit none
       real data(:)
       integer source
       integer, optional :: tag
       integer counter, given_tag

       given_tag=0
       if(present(tag)) given_tag=tag
       counter=size(data,1)
     
       call mpi_recv(data, counter*2, mpi_real, source, given_tag, MPI_COMM_ICE, mpistatus, ierr)

       if(ierr.ne.0) then
          print *, 'error common_mpi_recv_1D count=', counter, 'tag=', given_tag
          stop
       endif

      end subroutine common_mpi_recv_1D


       subroutine common_mpi_recv_2D(data, source, tag)
        use ice_communicate, only: MPI_COMM_ICE
        implicit none
        real data(:,:)
        integer source
        integer, optional :: tag
        integer counter, given_tag

        given_tag=0
        if(present(tag)) given_tag=tag
        counter=size(data,1)*size(data,2)

        call mpi_recv(data, counter, mpi_real, source, given_tag, MPI_COMM_ICE, mpistatus, ierr)

        if(ierr.ne.0) then
            print *, 'error common_mpi_recv_2D count=', counter, 'tag=', given_tag
            stop
        endif
       end subroutine common_mpi_recv_2D

end module common_mpi
