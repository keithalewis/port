      SUBROUTINE C6LCF(P,X,NF,F,IU,UR,UF)
      INTEGER P,IU
      REAL X(P),F,UR
      EXTERNAL UF
      CALL UF(P,X,NF,F)
      RETURN
      END
