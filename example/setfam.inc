  subroutine setfam(k, b)
    integer, intent(in) :: k
    type(fv_bicop), intent(out) :: b
    select case (k)
    case (1); b%family = FV_GAUSSIAN; b%par = [0.6_dp, 0.0_dp, 0.0_dp]
    case (2); b%family = FV_STUDENT; b%par = [0.5_dp, 4.5_dp, 0.0_dp]
    case (3); b%family = FV_CLAYTON; b%par = [2.5_dp, 0.0_dp, 0.0_dp]
    case (4); b%family = FV_GUMBEL; b%par = [2.2_dp, 0.0_dp, 0.0_dp]
    case (5); b%family = FV_FRANK; b%par = [6.0_dp, 0.0_dp, 0.0_dp]
    case (6); b%family = FV_JOE; b%par = [2.8_dp, 0.0_dp, 0.0_dp]
    case (7); b%family = FV_BB1; b%par = [0.8_dp, 1.6_dp, 0.0_dp]
    case (8); b%family = FV_BB6; b%par = [1.5_dp, 1.8_dp, 0.0_dp]
    case (9); b%family = FV_BB7; b%par = [1.6_dp, 1.2_dp, 0.0_dp]
    case (10); b%family = FV_BB8; b%par = [3.0_dp, 0.7_dp, 0.0_dp]
    case (11); b%family = FV_TAWN; b%par = [0.4_dp, 0.9_dp, 3.0_dp]
    case (12); b%family = FV_PLACKETT; b%par = [6.0_dp, 0.0_dp, 0.0_dp]
    case (13); b%family = FV_FGM; b%par = [0.7_dp, 0.0_dp, 0.0_dp]
    case (14); b%family = FV_AMH; b%par = [0.8_dp, 0.0_dp, 0.0_dp]
    case (15); b%family = FV_GALAMBOS; b%par = [1.4_dp, 0.0_dp, 0.0_dp]
    case (16); b%family = FV_HUSLER_REISS; b%par = [1.8_dp, 0.0_dp, 0.0_dp]
    case (17); b%family = FV_BB2; b%par = [0.8_dp, 1.2_dp, 0.0_dp]
    case (18); b%family = FV_BB3; b%par = [1.6_dp, 0.8_dp, 0.0_dp]
    case (19); b%family = FV_BB5; b%par = [1.8_dp, 0.9_dp, 0.0_dp]
    case (20); b%family = FV_BB9; b%par = [2.0_dp, 0.7_dp, 0.0_dp]
    case (21); b%family = FV_BB10; b%par = [2.0_dp, 0.8_dp, 0.0_dp]
    case (22); b%family = FV_BB4; b%par = [1.2_dp, 1.5_dp, 0.0_dp]
    case (23); b%family = FV_TEV; b%par = [0.5_dp, 4.0_dp, 0.0_dp]
    end select
  end subroutine setfam
