module M_testsuite_M_starpac
use M_verify
use M_starpac
character(len=*),parameter :: options=' -section 3 -library libGPF -filename `pwd`/M_starpac.FF &
& -documentation y -ufpp   y -ccall  n -archive  GPF.a '
contains
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
subroutine test_suite_M_starpac()
!   call test_adjustc()
!   call test_atleast()
!   call test_base()
end subroutine test_suite_M_starpac
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
subroutine test_base2()
character(len=:),allocatable :: in(:)
integer,allocatable          :: expected(:)
   call unit_check_start('base2',' '//OPTIONS)

   ! convert base10 values to base2 strings
   in=[character(len=32) :: '10','1010','101010','10101010','1010101010','101010101010']
   expected=[2,10,42,170,682,2730]
   call checkit(in,expected)

   call unit_check_done('base2')
contains
subroutine checkit(answer,values)
character(len=*),intent(in)  :: answer(:)
integer,intent(in)           :: values(:)
character(len=32)            :: out
integer                      :: i
   do i=1,size(answer)
      call unit_check('base2',base2(values(i)).eq.answer(i), &
       & 'checking for '//trim(answer(i))//' in base 2 from value '//v2s(values(i)) )
   enddo
end subroutine checkit
end subroutine test_base2
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
subroutine test_transliterate
!-!use M_starpac, only: transliterate
implicit none
character(len=36),parameter :: lc='abcdefghijklmnopqrstuvwxyz0123456789'
character(len=36),parameter :: uc='ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
   call unit_check_start('transliterate',' &
      & -description ''when characters in set one are found replace them with characters from set two'' '//OPTIONS )
!call unit_check('transliterate',transliterate('AbCDefgHiJklmnoPQRStUvwxyZ',lc,uc).eq.uc(1:26),msg='transliterate to uppercase')
!call unit_check('transliterate',transliterate('AbCDefgHiJklmnoPQRStUvwxyZ',uc,lc).eq.lc(1:26),msg='transliterate to lowercase')
call unit_check_done('transliterate')
end subroutine test_transliterate
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
end module M_testsuite_M_starpac
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
program runtest
use M_msg
use M_verify, only : unit_check_command, unit_check_keep_going, unit_check_level, unit_check_stop
use M_testsuite_M_starpac
   unit_check_command=''
   unit_check_keep_going=.true.
   unit_check_level=0
   call test_suite_M_starpac()
   call unit_check_stop()
end program runtest
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
