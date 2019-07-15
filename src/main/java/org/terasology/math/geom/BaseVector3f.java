/*
 * Copyright 2019 MovingBlocks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.terasology.math.geom;

import com.google.common.base.Preconditions;
import org.joml.Matrix3dc;
import org.joml.Matrix3fc;
import org.joml.Matrix3x2fc;
import org.joml.Matrix4dc;
import org.joml.Matrix4fc;
import org.joml.Matrix4x3fc;
import org.joml.Quaternionf;
import org.joml.Quaternionfc;
import org.joml.Vector3fc;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;

/**
 * A vector/point in 3D space
 */
public abstract class BaseVector3f implements Vector3fc{

    /**
     * An immutable instance with all components set to 0
     */
    public static final ImmutableVector3f ZERO = new ImmutableVector3f(0, 0, 0);

    /**
     * An immutable instance with all components set to 1
     */
    public static final ImmutableVector3f ONE = new ImmutableVector3f(1, 1, 1);

    /**
      * @return x the x coordinate
      */
    public abstract float getX(); 
    /**
      * @return y the y coordinate
      */
    public abstract float getY(); 
    /**
      * @return z the z coordinate
      */
    public abstract float getZ(); 

    /**
      * @return x the x coordinate
      */
    public abstract float x(); 
    /**
      * @return y the y coordinate
      */
    public abstract float y(); 
    /**
      * @return z the z coordinate
      */
    public abstract float z(); 

    /**
     * @param a the first point
     * @param b the second point
     * @param t the interpolation value in the range [0..1]
     * @return the interpolated point
     */
    public static Vector3f lerp(Vector3fc a, Vector3fc b, float t) {
        Preconditions.checkArgument(t >= 0 && t <= 1, "t must be in range [0..1]");

        float x = a.x() * (1 - t) + b.x() * t;
        float y = a.y() * (1 - t) + b.y() * t;
        float z = a.z() * (1 - t) + b.z() * t;
        return new Vector3f(x, y, z);
    }
    /**
     * Returns the dot product of this vector and vector other.
     * @param other the other vector
     * @return the dot product of this and other
     */
//    public final float dot(Vector3fc other) {
//        return (float) (this.x() * other.x() + this.y() * other.y() + this.z() * other.z());
//    }


    /**
    * Returns the current vector projected onto v
    * @param v other vector to project onto
    *
    */
    public final Vector3f project( Vector3fc v)
    {
        return new Vector3f(v).mul(this.dot(v)/ (v.lengthSquared()));
    }

    /**
     * @return the squared distance to the origin
     */
    public float lengthSquared() {
        return getX() * getX() + getY() * getY() + getZ() * getZ();
    }
   /**
    *   Returns the angle in radians between this vector and the vector
    *   parameter; the return value is constrained to the range [0,PI].
    *   @param v1    the other vector
    *   @return   the angle in radians in the range [0,PI]
    */
//   public final float angle(Vector3fc v1) {
//      double vDot = this.dot(v1) / (this.length() * v1.length());
//
//      if (vDot < -1.0) {
//          vDot = -1.0;
//      }
//
//      if (vDot >  1.0) {
//          vDot =  1.0;
//      }
//
//      return (float) Math.acos(vDot);
//   }

    /**
     * @return the distance to the origin
     */
    public float length() {
        return (float) Math.sqrt(lengthSquared());
    }

    /**
     * @param other the other point
     * @return the distance in between
     */
//    public float distanceSquared(Vector3fc other) {
//        float dx = other.x() - this.getX();
//        float dy = other.y() - this.getY();
//        float dz = other.z() - this.getZ();
//
//        return dx * dx + dy * dy + dz * dz;
//    }

    /**
     * @param other the other point
     * @return the distance in between
     */
//    public float distance(Vector3fc other) {
//        return (float) Math.sqrt(distanceSquared(other));
//    }


    /**
     * Computes the distance between two points
     * @param p1 the first point
     * @param p2 the second point
     * @return the distance between the two points
     */
    public static float distance(BaseVector3f p1, BaseVector3f p2) {
        return p1.distance(p2);
    }

    /**
     * All point implementations with the same coordinate are equal
     */
    @Override
    public final boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj instanceof BaseVector3f) {
            BaseVector3f other = (BaseVector3f) obj;
            return Float.floatToIntBits(getX()) == Float.floatToIntBits(other.getX())
                && Float.floatToIntBits(getY()) == Float.floatToIntBits(other.getY())
                && Float.floatToIntBits(getZ()) == Float.floatToIntBits(other.getZ());
        }
        return false;
    }

    /**
     * All point implementations with the same coordinate have the same hashcode
     */
    @Override
    public final int hashCode() {
        int result = 1;
        final int prime = 31;
        long temp;
        temp = Float.floatToIntBits(getX());
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Float.floatToIntBits(getY());
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Float.floatToIntBits(getZ());
        result = prime * result + (int) (temp ^ (temp >>> 32));

        return result;
    }

    @Override
    public String toString() {
        return "(" + getX() + ", " + getY() + ", " + getZ() + ")";
    }

    @Override
    public ByteBuffer get(ByteBuffer buffer) {
        return new Vector3f(this).get(buffer);
    }

    @Override
    public ByteBuffer get(int index, ByteBuffer buffer) {
        return new Vector3f(this).get(index,buffer);
    }

    @Override
    public FloatBuffer get(FloatBuffer buffer) {
        return new Vector3f(this).get(buffer);
    }

    @Override
    public FloatBuffer get(int index, FloatBuffer buffer) {
        return new Vector3f(this).get(index,buffer);
    }

    @Override
    public Vector3fc  getToAddress(long address) {
        return new Vector3f(this).getToAddress(address);
    }

    @Override
    public org.joml.Vector3f mul(float scalar, org.joml.Vector3f dest) {
       dest.x = scalar * getX();
       dest.y = scalar * getY();
       dest.z = scalar * getZ();
       return dest;
    }

    @Override
    public org.joml.Vector3f mul(Vector3fc v, org.joml.Vector3f dest) {
        dest.x = v.x() * x();
        dest.y = v.y() * y();
        dest.z = v.z() * z();

        return dest;
    }

    @Override
    public org.joml.Vector3f div(Vector3fc v, org.joml.Vector3f dest) {
        return new Vector3f(this).div(v,dest);
    }

    @Override
    public float distanceSquared(Vector3fc v) {
        return new org.joml.Vector3f(this).distanceSquared(v);
    }

    @Override
    public float distanceSquared(float x, float y, float z) {
       return new org.joml.Vector3f(this).distanceSquared(x,y,z);
    }
    @Override
    public float dot(Vector3fc v) {
        return x() * v.x() +y() * v.y() +z() * v.z() ;
    }

    @Override
    public float dot(float x, float y, float z) {
        return x()*x + y()*y + z()*z;
    }

    @Override
    public float angleCos(Vector3fc v) {
        return new Vector3f(this).angleCos(v);
    }

    @Override
    public float angle(Vector3fc v) {
        return new org.joml.Vector3f(this).angle(v);
    }

    @Override
    public org.joml.Vector3f normalize(org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).normalize(dest);
    }

    @Override
    public org.joml.Vector3f normalize(float length, org.joml.Vector3f dest) {
         return new org.joml.Vector3f(this).normalize(length,dest);
    }

    @Override
    public org.joml.Vector3f mul(float x, float y, float z, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(x,y,z).mul(dest);
    }


    @Override
    public  float distance(Vector3fc v) {
        return new org.joml.Vector3f(this).distance(v);
    }

    @Override
    public float distance(float x, float y, float z) {
       return new org.joml.Vector3f(this).distance(x,y,z);
    }


    @Override
    public org.joml.Vector3f sub(Vector3fc v, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).sub(v,dest);
    }


    @Override
    public org.joml.Vector3f sub(float x, float y, float z, org.joml.Vector3f dest) {
        dest.x = x() - x;
        dest.y = y() - y;
        dest.z = z() - z;
        return dest;
    }


    @Override
    public org.joml.Vector3f add(Vector3fc v, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).add(v,dest);
    }

    @Override
    public  org.joml.Vector3f add(float x, float y, float z,  org.joml.Vector3f dest) {
        dest.x = x() + x;
        dest.y = y() + y;
        dest.z = z() + z;
        return dest;
    }

    @Override
    public org.joml.Vector3f negate(org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).negate(dest);
    }

    @Override
    public org.joml.Vector3f absolute(org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).absolute(dest);
    }

    @Override
    public org.joml.Vector3f reflect(Vector3fc normal, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).reflect(normal,dest);
    }

    @Override
    public org.joml.Vector3f reflect(float x, float y, float z, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).reflect(x,y,z,dest);
    }

    @Override
    public org.joml.Vector3f half(Vector3fc other, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).half(other,dest);
    }

    @Override
    public org.joml.Vector3f half(float x, float y, float z, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).half(x,y,z,dest);
    }

    @Override
    public org.joml.Vector3f lerp(Vector3fc other, float t, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).lerp(other,t,dest);
    }
    @Override
    public org.joml.Vector3f fma(Vector3fc a, Vector3fc b, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).fma(b,dest);
    }

    @Override
    public org.joml.Vector3f fma(float a, Vector3fc b, org.joml.Vector3f dest) {
      return new org.joml.Vector3f(this).fma(a,b,dest);
    }

    @Override
    public org.joml.Vector3f min(Vector3fc v, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).min(v,dest);
    }

    @Override
    public org.joml.Vector3f max(Vector3fc v, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).max(v,dest);
    }

    @Override
    public float get(int component) throws IllegalArgumentException {
       return new org.joml.Vector3f(this).get(component);
    }

    @Override
    public org.joml.Vector3f mulProject(Matrix4fc mat, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).mulProject(mat,dest);
    }

    @Override
    public org.joml.Vector3f mul(Matrix3dc mat, org.joml.Vector3f dest) {
      return new org.joml.Vector3f(this).mul(mat,dest);
    }

    @Override
    public org.joml.Vector3f mul(Matrix3fc mat, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).mul(mat,dest);
    }
    @Override
    public org.joml.Vector3f mul(Matrix3x2fc mat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).mul(mat,dest);
    }
    @Override
    public org.joml.Vector3f mulTranspose(Matrix3fc mat, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).mulTranspose(mat,dest);
    }
    @Override
    public org.joml.Vector3f mulPosition(Matrix4fc mat, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).mulPosition(mat,dest);
    }
    @Override
    public org.joml.Vector3f mulPosition(Matrix4x3fc mat, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).mulPosition(mat,dest);
    }

    @Override
    public org.joml.Vector3f mulTransposePosition(Matrix4fc mat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).mulTransposePosition(mat,dest);
    }

    @Override
    public float mulPositionW(Matrix4fc mat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).mulPositionW(mat,dest);
    }

    @Override
    public org.joml.Vector3f mulDirection(Matrix4dc mat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).mulDirection(mat,dest);
    }

    @Override
    public org.joml.Vector3f mulDirection(Matrix4fc mat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).mulDirection(mat,dest);
    }

    @Override
    public org.joml.Vector3f mulDirection(Matrix4x3fc mat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).mulDirection(mat,dest);
    }

    @Override
    public org.joml.Vector3f mulTransposeDirection(Matrix4fc mat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).mulTransposeDirection(mat,dest);
    }

    @Override
    public org.joml.Vector3f rotate(Quaternionfc quat, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).rotate(quat,dest);
    }

    @Override
    public Quaternionf rotationTo(Vector3fc toDir, Quaternionf dest) {
        return new org.joml.Vector3f(this).rotationTo(toDir,dest);
    }

    @Override
    public Quaternionf rotationTo(float toDirX, float toDirY, float toDirZ, Quaternionf dest) {
        return new org.joml.Vector3f(this).rotationTo(toDirX,toDirY,toDirZ,dest);
    }

    @Override
    public org.joml.Vector3f rotateAxis(float angle, float aX, float aY, float aZ, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).rotateAxis(angle,aX,aY,aZ,dest);
    }

    @Override
    public org.joml.Vector3f rotateX(float angle, org.joml.Vector3f dest) {
         return new org.joml.Vector3f(this).rotateX(angle,dest);
    }

    @Override
    public org.joml.Vector3f rotateY(float angle, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).rotateY(angle,dest);
    }

    @Override
    public org.joml.Vector3f rotateZ(float angle, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).rotateZ(angle,dest);
    }

    @Override
    public org.joml.Vector3f div(float scalar, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).div(scalar,dest);
    }

    @Override
    public org.joml.Vector3f div(float x, float y, float z, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).div(x,y,z,dest);
    }

    @Override
    public org.joml.Vector3f cross(Vector3fc v, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).cross(v,dest);
    }

    @Override
    public org.joml.Vector3f cross(float x, float y, float z, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).cross(x,y,z,dest);
    }

    @Override
    public org.joml.Vector3f smoothStep(Vector3fc v, float t, org.joml.Vector3f dest) {
       return new org.joml.Vector3f(this).smoothStep(v,t,dest);
    }

    @Override
    public org.joml.Vector3f hermite(Vector3fc t0, Vector3fc v1, Vector3fc t1, float t, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).hermite(t0,v1,t1,t,dest);
    }

    @Override
    public int maxComponent() {
        return new org.joml.Vector3f(this).maxComponent();
    }

    @Override
    public int minComponent() {
        return new org.joml.Vector3f(this).minComponent();
    }

    @Override
    public org.joml.Vector3f orthogonalize(Vector3fc v, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).orthogonalize(v,dest);
    }

    @Override
    public org.joml.Vector3f orthogonalizeUnit(Vector3fc v, org.joml.Vector3f dest) {
        return new org.joml.Vector3f(this).orthogonalizeUnit(v,dest);
    }


    @Override
    public boolean equals(Vector3fc v, float delta) {
        return new org.joml.Vector3f(v.x(),v.y(),v.z()).equals(v,delta);
    }


}
