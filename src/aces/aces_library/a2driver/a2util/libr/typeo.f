
C CHARACTER*8 FUNCTION WHICH RETURNS THE TYPE OF AN ORBITAL (OCCUPIED
C OR VIRTUAL) FOR AN GIVEN ORBITAL BY CHECKING ITS INDEX

      CHARACTER*8 FUNCTION TYPEO(I,NI)
      IF (I.LE.NI) THEN
         TYPEO='occupied'
      ELSE
         TYPEO=' virtual'
      END IF
      RETURN
      END