!     this program was adapted with help from 'pmy' on stack overflow
!     URL: https://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90
!
!
!
!
      program memory
!
      implicit none
!
      integer :: valueRSS
      character(len=200):: filename=' '
      character(len=80) :: line
      character(len=8)  :: pid_char=' '
      integer :: pid
      logical :: ifxst
!
!
      integer :: i, tmax, sle
      real :: time0, timei
!
!
      valueRSS=-1
      pid = 9999999
      tmax = 1200
      sle = 1
!
!
!
!
!
!
!
!     remove traces of old programs
      call EXECUTE_COMMAND_LINE('rm PROCESS')
!
!
!
!
!
!
!
!     check to see if other program has started
      do i = 1,tmax
!
      call SLEEP(1)
!
      inquire (file='PROCESS',exist=ifxst)
      if (.not.ifxst) then
        print*, "Scanning for 'PROCESS' file . . ."
      else
        print*, "Found file."
        call CPU_TIME(time0)
        go to 123
      endif
!
      if (i .EQ. tmax) then
        print*, "Scanning timed out."
        return
      end if
!
      end do
123   continue
!
!
!
!
!
!
!     if it has, read the 'process' file to get the process ID
!
      open(unit=200, file='PROCESS', action='read')
      read (200,*) pid
      close(200)
!
      print*, "Found PID to use:", pid
      write(pid_char,'(I8)') pid
      filename='/proc/'//trim(adjustl(pid_char))//'/status'
!
!
!
!
!
!     open file to record results
      open(888, file='memory.plt')
      write(888,*) "Initial time", time0
      write(888,*) 'VARIABLES = "INDEX" "TIME" "MEMORY"'
!
!
!
!
!
!
!
!
!
!
!     scan the file until escape conditions
      do i = 1,tmax
!
!
      call SLEEP(1)
!
!
      inquire (file=filename,exist=ifxst)
      if (.not.ifxst) then
        write (*,*) 'Monitored program ended. Exiting program.'
        close(888)
        return
      endif
!
!
      call CPU_TIME(timei)
      open(unit=100, file=filename, action='read')
      do
        read (100,'(a)',end=120) line
        if (line(1:6).eq.'VmRSS:') then
           read (line(7:),*) valueRSS
           exit
        endif
      enddo
120   continue
      close(100)
!
      print*, i, timei-time0, "RSS value:", valueRSS
      write(888,*) i, timei-time0, valueRSS
!
!
!
      end do
!
!
!
!
!
!
!
!
!
!
!
!
      end program memory 
