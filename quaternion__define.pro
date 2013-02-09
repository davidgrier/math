;+
; NAME:
;    Quaternion__define
;
; PURPOSE:
;    Class for quaternions
;
; CATEGORY:
;    Mathematics, objects
;
; PROPERTIES:
;    X, Y, Z: Imaginary (vector) components of quaternion
;    W:       Real (scalar) component 
;
;    In describing rotations, (x,y,z) corresponds to the axis of
;    rotation and w is the rotation angle.
;
; METHODS:
;    Quaternion::getproperty
;
;    Quaternion::setproperty
;
; OVERLOADED OPERATORS:
;    ~: Quaternion conjugate.  Negates imaginary part.
;    +: Component-by-component addition
;    -: Component-by-component subtraction
;    *: Multiplication by scalars, vectors, quaternions,
;        or arrays of these.  NOTE: not all implemented yet
;    /: Division by scalars, vectors, quaternions,
;        or arrays of these.  NOTE: not all implemented yet
;
; MODIFICATION HISTORY:
; 01/29/2011 Written by David G. Grier, New York University
; 02/04/2011 DGG overloaded subscripting to access individual components
;
; Copyright (c) 2011, David G. Grier
;-    

;;;;;
;
; quaternion::rotatevector(v)
;
; Rotate vector v using a quaternion
;
function quaternion::RotateVector, v

if n_elements(v) eq 3 then begin
   res = self * v * (~self)
   res.getproperty, x = x, y = y, z = z
   return, [x, y, z]
endif

return, v
end

;;;;;
;
; OPERATOR OVERLOADING
;
;;;;;
;
; Plus
;
function quaternion::_overloadPlus, left, right

left->getproperty, x = x1, y = y1, z = z1, w = w1
right->getproperty, x = x2, y = y2, z = z2, w = w2
return, quaternion(x = x1+x2, y = y1+y2, z = z1+z2, w = w1+w2)
end

;;;;;
;
; Minus
;
function quaternion::_overloadMinus, left, right

left->getproperty, x = x1, y = y1, z = z1, w = w1
right->getproperty, x = x2, y = y2, z = z2, w = w2
return, quaternion(x = x1-x2, y = y1-y2, z = z1-z2, w = w1-w2)
end

;;;;;
;
; Multiply
;
function quaternion::_overloadAsterisk, left, right

if isa(left, "quaternion") then begin
   q = left.duplicate()
   q.multiply, right
endif else if n_elements(left) eq 3 then begin
   q = quaternion(vector = float(left), angle = 0.)
   q.multiply, right
endif else if n_elements(left) eq 1 then begin
   q = right.duplicate()
   q.multiply, left
endif

return, q
end

;;;;;
;
; Divide
;
function quaternion::_overloadSlash, q1, q2

q1->getproperty, x = x1, y = y1, z = z1, w = w1
q2->getproperty, x = x2, y = y2, z = z2, w = w2
norm = q2.Norm()
w = (w2 * w1 + x2 * x1 + y2 * y1 + z2 * z1) / norm
x = (w2 * x1 - x2 * w1 - y2 * z1 + z2 * y1) / norm
y = (w2 * y1 + x2 * z1 - y2 * w1 - z2 * x1) / norm
z = (w2 * z1 - x2 * y1 + y2 * x1 - z2 * w1) / norm
return, quaternion(w = w, x = x, y = y, z = z)
end

;;;;;
;
; Conjugate (Tilde)
;
function quaternion::_overloadTilde

return, quaternion(w = self.w, x = -self.x, y = -self.y, z = -self.z)
end

;;;;;
;
; Subscripting
;
function quaternion::_overloadBracketsRightSide, isrange, sub1

on_error, 2
if isrange eq 0 then begin
   case sub1 of
      0: return, self.x
      1: return, self.y
      2: return, self.z
      3: return, self.w
      else: message, 'subscript out of range'
   endcase
endif

temp = [self.x, self.y, self.z, self.w]
return, temp[sub1[0]:sub1[1]:sub1[2]]
end

pro quaternion::_overloadBracketsLeftSide, objref, rightvalue, isrange, sub1

on_error, 2

if isrange eq 0 then begin
   case sub1 of
      0: self.x = rightvalue
      1: self.y = rightvalue
      2: self.z = rightvalue
      3: self.w = rightvalue
      else: message, 'subscript out of range'
   endcase
endif else begin
   count = 0
   for i = sub1[0], sub1[1], sub1[2] do begin
      case i of
         0: self.x = rightvalue[count++]
         1: self.y = rightvalue[count++]
         2: self.z = rightvalue[count++]
         3: self.w = rightvalue[count++]
         else: message, 'subscript out of range'
      endcase
   endfor
endelse

end

;;;;;
;
; Print
;
function quaternion::_overloadPrint

return, [self.x, self.y, self.z, self.w]
end

;;;;;
;
; Help
;
function quaternion::_overloadHelp, varname

type = "QUATERNION"
value = "("+strjoin(strtrim([self.x, self.y, self.z, self.w], 2), ', ')+")"
format = '(%"%-15s %-10s = %s")'

return, string(varname, type, value, format = format)
end

;;;;;
;
; quaternion::Multiply
;
; result = self * this: left multiplication
;
pro quaternion::Multiply, this

if isa(this, 'quaternion') then $
   this->getproperty, x = x, y = y, z = z, w = w $
else if n_elements(this) eq 3 then begin
   x = this[0]
   y = this[1]
   z = this[2]
   w = 0.
endif else if n_elements(this) eq 1 then begin
   self.x *= this
   self.y *= this
   self.z *= this
   self.w *= this
   return
endif
;q1->getproperty, x = x1, y = y1, z = z1, w = w1
;q2->getproperty, x = x2, y = y2, z = z2, w = w2
;w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
;x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
;y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
;z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2

tw = self.w * w - self.x * x - self.y * y - self.z * z
tx = self.w * x + self.x * w + self.y * z - self.z * y
ty = self.w * y - self.x * z + self.y * w + self.z * x
tz = self.w * z + self.x * y - self.y * x + self.z * w

self.w = tw
self.x = tx
self.y = ty
self.z = tz
end

;;;;;
;
; quaternion::Divide
;
; Divide one quaternion by another
;
pro quaternion::Divide, this

if ~isa(this, 'quaternion') then begin
   message, 'argument must be a quaternion', /inf
   return
endif

this->getproperty, x = x, y = y, z = z, w = w
norm = this.Norm()
tw = w * self.w + x * self.x + y * self.y + z * self.z
tx = w * self.x - x * self.w - y * self.z + z * self.y
ty = w * self.y + x * self.z - y * self.w - z * self.x
tz = w * self.z - x * self.y + y * self.x - z * self.w
self.w = tw / norm
self.x = tx / norm
self.y = ty / norm
self.z = tz / norm
end

;;;;;
;
; quaternion::Add
;
; Add to a quaternion
;
pro quaternion::Add, this

; add scalar quaternion
if isa(this, 'quaternion') then begin
   this.getproperty, x = dx, y = dy, z = dz, w = dw
   self.x += dx
   self.y += dy
   self.z += dz
   self.w += dw
endif

end

;;;;;
;
; quaternion::Inverse
;
; Compute inverse of quaternion
;
pro quaternion::Inverse

self.conjugate
len = self.norm()
self.x /= len
self.y /= len
self.z /= len
self.w /= len
end

;;;;;
;
; quaternion::Conjugate
;
; Take the conjugate of a quaternion
;
pro quaternion::Conjugate

self.x *= -1.
self.y *= -1.
self.z *= -1.
end

;;;;;
;
; quaternion::Modulus()
;
; Compute modulus of quaternion
;
function quaternion::Modulus

return, sqrt(self.norm())
end

;;;;;
;
; quaternion::Norm()
;
; Norm of a quaternion
;
function quaternion::Norm

return, self.x^2 + self.y^2 + self.z^2 + self.w^2
end

;;;;;
;
; quaternion::Normalize
;
; Normalize a quaternion
;
pro quaternion::Normalize

len = self.modulus()

if len lt 1e-5 then return

self.x /= len
self.y /= len
self.z /= len
self.w /= len
end

;;;;;
;
; quaternion::Duplicate
;
; Create a dupblicate of the present quaternion
;
function quaternion::Duplicate

return, quaternion(x = self.x, y = self.y, z = self.z, w = self.w)
end

;;;;;
;
; quaternion::GetProperty
;
; Get quaternion properties
;
pro quaternion::GetProperty, x = x, $
                             y = y, $
                             z = z, $
                             w = w

x = self.x
y = self.y
z = self.z
w = self.w
end

;;;;;
;
; quaternion::SetProperty
;
; Set quaternion properties
;
pro quaternion::SetProperty, x = x, $
                             y = y, $
                             z = z, $
                             w = w, $
                             vector = vector, $
                             angle = angle, $
                             v1 = v1, $
                             v2 = v2

if n_elements(x) eq 1 then self.x = float(x)
if n_elements(y) eq 1 then self.y = float(y)
if n_elements(z) eq 1 then self.z = float(z)
if n_elements(w) eq 1 then self.w = float(w)

if n_elements(vector) eq 3 then begin
   self.x = float(vector[0])
   self.y = float(vector[1])
   self.z = float(vector[2])
endif

if n_elements(w) eq 1 then self.w = float(w)

if n_elements(angle) eq 1 then self.w = float(angle)

if (n_elements(v1) eq 3) and (n_elements(v2) eq 3) then begin
   v = crossp(v1, v2)
   self.x = v[0]
   self.y = v[1]
   self.z = v[2]
   self.w = sqrt(total(v1^2)) * sqrt(total(v2^2)) + total(v1 * v2)
endif

end

;;;;;
;
; quaternion::Cleanup
;
; Free resources
;
pro quaternion::Cleanup

end

;;;;;
;
; quaternion::Init
;
; Initialize the quaternion object either as a vector and an angle,
; or else with a starting vector (v1) and an ending vector (v2)
;
function quaternion::Init, x = x,           $
                           y = y,           $
                           z = z,           $
                           w = w,           $
                           vector = vector, $
                           angle = angle,   $
                           v1 = v1,         $
                           v2 = v2

v = fltarr(3)

if n_elements(x) eq 1 then v[0] = float(x)
if n_elements(y) eq 1 then v[1] = float(y)
if n_elements(z) eq 1 then v[2] = float(z)
if n_elements(w) eq 1 then self.w = float(w)

if n_elements(vector) eq 3 then $
   v = float(vector)

if n_elements(angle) eq 1 then $
   self.w = float(angle)

if (n_elements(v1) eq 3) and (n_elements(v2) eq 3) then begin
   v = crossp(v1, v2)
   self.w = sqrt(total(v1^2)) * sqrt(total(v2^2)) + total(v1 * v2)
endif

self.x = v[0]
self.y = v[1]
self.z = v[2]

return, 1
end

;;;;;
;
; quaternion__define
;
; Define the quaternion object
;
pro quaternion__define

struct = {quaternion, $
          inherits IDL_Object, $
          x: 0., $
          y: 0., $
          z: 0., $
          w: 0.  $
         }
end
