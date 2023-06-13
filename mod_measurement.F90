module mod_measurement

   use ice_kinds_mod

   type measurement
      real d                       ! Measurement value
      real var                     ! Error variance of measurement
      character(len=3) id          ! Type of measurement (ex)SIT)
      real lon                     ! Longitude position
      real lat                     ! Latitude position
      integer ipiv                 ! i-pivot point in grid
      integer jpiv                 ! j-pivot point in grid
      real a1                      ! bilinear coeffisients (if ni=0)
      real a2                      ! bilinear coeffisients
      real a3                      ! bilinear coeffisients
      real a4                      ! bilinear coeffisients

   end type measurement

   type inimeasurement4
      real*4 d                       ! Measurement value
      real*4 var                     ! Error variance of measurement
      character(len=3) id            ! Type of measurement (ex)SIT)
      real*4 lon                     ! Longitude position
      real*4 lat                     ! Latitude position

   end type inimeasurement4

end module mod_measurement

