!// ForQuill v1.01 Beta www.fcode.cn
Subroutine ptek(e, l, t, np, ne)
  Dimension x(50), y(50), a(80), nenp(80, 2), ek(4, 4), t(4, 4), s(4, 4)
  Common nenp, x, y, a, ek
  i = nenp(l, 1)
  j = nenp(l, 2)
  xi = x(i)
  xj = x(j)
  yi = y(i)
  yj = y(j)
  al = sqrt((xj-xi)**2+(yj-yi)**2)
  Do i = 1, 4
    Do j = 1, 4
      ek(i, j) = 0.
    End Do
  End Do
  c = a(l)*e/al
  ek(1, 1) = c
  ek(1, 3) = -c
  ek(3, 1) = -c
  ek(3, 3) = c
  Do i = 1, 4
    Do j = 1, 4
      t(i, j) = 0.
    End Do
  End Do
  sn = (yj-yi)/al
  cs = (xj-xi)/al
  t(1, 1) = cs
  t(1, 2) = -sn
  t(2, 1) = sn
  t(2, 2) = cs
  Do i = 1, 2
    Do j = 1, 2
      t(i+2, j+2) = t(i, j)
    End Do
  End Do
  Do i = 1, 4
    Do j = 1, 4
      s(i, j) = 0.
      Do k = 1, 4
        s(i, j) = s(i, j) + t(i, k)*ek(k, j)
      End Do
    End Do
  End Do
  Do i = 1, 4
    Do j = 1, 4
      ek(i, j) = 0.
      Do k = 1, 4
        ek(i, j) = ek(i, j) + s(i, k)*t(j, k)
      End Do
    End Do
  End Do
  Return
End Subroutine ptek