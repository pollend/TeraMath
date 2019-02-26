/*
 * Copyright 2019 MovingBlocks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.terasology.math.geom;

import org.joml.AxisAngle4d;
import org.joml.AxisAngle4f;
import org.joml.Matrix3d;
import org.joml.Matrix3dc;
import org.joml.Matrix3fc;
import org.joml.Quaterniond;
import org.joml.Quaterniondc;
import org.joml.Quaternionf;
import org.joml.Quaternionfc;
import org.joml.Vector3d;
import org.joml.Vector3dc;
import org.joml.Vector3f;
import org.joml.Vector3fc;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;

/**
 * Defines an immutable 3x3 double matrix
 * @author auto-generated
 */
public class ImmutableMatrix3d extends BaseMatrix3d {

    org.joml.Matrix3d matrix;

    /**
     * Constructs and initializes a Matrix3d from the specified values.
     * @param m00 the m00 component
     * @param m01 the m01 component
     * @param m02 the m02 component
     * @param m10 the m10 component
     * @param m11 the m11 component
     * @param m12 the m12 component
     * @param m20 the m20 component
     * @param m21 the m21 component
     * @param m22 the m22 component
     */
    public ImmutableMatrix3d(double m00, double m01, double m02, double m10, double m11, double m12, double m20, double m21, double m22) {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
    }

    /**
     *  Constructs a new matrix with the same values as the
     *  Matrix3d parameter.
     *  @param m1  the source matrix
     */
    public ImmutableMatrix3d(BaseMatrix3d m1) {
        this.m00 = m1.getM00();
        this.m01 = m1.getM01();
        this.m02 = m1.getM02();
        this.m10 = m1.getM10();
        this.m11 = m1.getM11();
        this.m12 = m1.getM12();
        this.m20 = m1.getM20();
        this.m21 = m1.getM21();
        this.m22 = m1.getM22();
    }

    @Override
    public final double getM00() {
        return m00;
    }

    @Override
    public final double getM01() {
        return m01;
    }

    @Override
    public final double getM02() {
        return m02;
    }

    @Override
    public final double getM10() {
        return m10;
    }

    @Override
    public final double getM11() {
        return m11;
    }

    @Override
    public final double getM12() {
        return m12;
    }

    @Override
    public final double getM20() {
        return m20;
    }

    @Override
    public final double getM21() {
        return m21;
    }

    @Override
    public final double getM22() {
        return m22;
    }


    @Override
    public double m00() {
        return 0;
    }

    @Override
    public double m01() {
        return 0;
    }

    @Override
    public double m02() {
        return 0;
    }

    @Override
    public double m10() {
        return 0;
    }

    @Override
    public double m11() {
        return 0;
    }

    @Override
    public double m12() {
        return 0;
    }

    @Override
    public double m20() {
        return 0;
    }

    @Override
    public double m21() {
        return 0;
    }

    @Override
    public double m22() {
        return 0;
    }

    @Override
    public Matrix3d mul(Matrix3dc right, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d mulLocal(Matrix3dc left, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d mul(Matrix3fc right, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d invert(Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d transpose(Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d get(Matrix3d dest) {
        return null;
    }

    @Override
    public AxisAngle4f getRotation(AxisAngle4f dest) {
        return null;
    }

    @Override
    public Quaternionf getUnnormalizedRotation(Quaternionf dest) {
        return null;
    }

    @Override
    public Quaternionf getNormalizedRotation(Quaternionf dest) {
        return null;
    }

    @Override
    public Quaterniond getUnnormalizedRotation(Quaterniond dest) {
        return null;
    }

    @Override
    public Quaterniond getNormalizedRotation(Quaterniond dest) {
        return null;
    }

    @Override
    public DoubleBuffer get(DoubleBuffer buffer) {
        return null;
    }

    @Override
    public DoubleBuffer get(int index, DoubleBuffer buffer) {
        return null;
    }

    @Override
    public FloatBuffer get(FloatBuffer buffer) {
        return null;
    }

    @Override
    public FloatBuffer get(int index, FloatBuffer buffer) {
        return null;
    }

    @Override
    public ByteBuffer get(ByteBuffer buffer) {
        return null;
    }

    @Override
    public ByteBuffer get(int index, ByteBuffer buffer) {
        return null;
    }

    @Override
    public ByteBuffer getFloats(ByteBuffer buffer) {
        return null;
    }

    @Override
    public ByteBuffer getFloats(int index, ByteBuffer buffer) {
        return null;
    }

    @Override
    public Matrix3dc getToAddress(long address) {
        return null;
    }

    @Override
    public double[] get(double[] arr, int offset) {
        return new double[0];
    }

    @Override
    public float[] get(float[] arr, int offset) {
        return new float[0];
    }

    @Override
    public float[] get(float[] arr) {
        return new float[0];
    }

    @Override
    public Matrix3d scale(Vector3dc xyz, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d scale(double x, double y, double z, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d scale(double xyz, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d scaleLocal(double x, double y, double z, Matrix3d dest) {
        return null;
    }

    @Override
    public Vector3d transform(Vector3d v) {
        return null;
    }

    @Override
    public Vector3d transform(Vector3dc v, Vector3d dest) {
        return null;
    }

    @Override
    public Vector3f transform(Vector3f v) {
        return null;
    }

    @Override
    public Vector3f transform(Vector3fc v, Vector3f dest) {
        return null;
    }

    @Override
    public Vector3d transform(double x, double y, double z, Vector3d dest) {
        return null;
    }

    @Override
    public Vector3d transformTranspose(Vector3d v) {
        return null;
    }

    @Override
    public Vector3d transformTranspose(Vector3dc v, Vector3d dest) {
        return null;
    }

    @Override
    public Vector3d transformTranspose(double x, double y, double z, Vector3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateX(double ang, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateY(double ang, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateZ(double ang, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateXYZ(double angleX, double angleY, double angleZ, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateZYX(double angleZ, double angleY, double angleX, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateYXZ(double angleY, double angleX, double angleZ, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotate(double ang, double x, double y, double z, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateLocal(double ang, double x, double y, double z, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateLocalX(double ang, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateLocalY(double ang, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateLocalZ(double ang, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateLocal(Quaterniondc quat, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateLocal(Quaternionfc quat, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotate(Quaterniondc quat, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotate(Quaternionfc quat, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotate(AxisAngle4f axisAngle, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotate(AxisAngle4d axisAngle, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotate(double angle, Vector3dc axis, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotate(double angle, Vector3fc axis, Matrix3d dest) {
        return null;
    }

    @Override
    public Vector3d getRow(int row, Vector3d dest) throws IndexOutOfBoundsException {
        return null;
    }

    @Override
    public Vector3d getColumn(int column, Vector3d dest) throws IndexOutOfBoundsException {
        return null;
    }

    /**
     * Retrieves the value at the specified row and column of the specified
     * matrix.
     * @param row the row number to be retrieved (zero indexed)
     * @param column the column number to be retrieved (zero indexed)
     * @return the value at the indexed element.
     */
    @Override
    public final double get(int row, int column) {
        switch (row) {
            case 0:
                switch (column) {
                    case 0:
                        return (this.getM00());
                    case 1:
                        return (this.getM01());
                    case 2:
                        return (this.getM02());
                    default:
                        break;
                }
                break;
            case 1:
                switch (column) {
                    case 0:
                        return (this.getM10());
                    case 1:
                        return (this.getM11());
                    case 2:
                        return (this.getM12());
                    default:
                        break;
                }
                break;

            case 2:
                switch (column) {
                    case 0:
                        return (this.getM20());
                    case 1:
                        return (this.getM21());
                    case 2:
                        return (this.getM22());
                    default:
                        break;
                }
                break;

            default:
                break;
        }

        throw new ArrayIndexOutOfBoundsException("row/col not in [0..2]");
    }

    @Override
    public Matrix3d normal(Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d lookAlong(Vector3dc dir, Vector3dc up, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d lookAlong(double dirX, double dirY, double dirZ, double upX, double upY, double upZ, Matrix3d dest) {
        return null;
    }

    @Override
    public Vector3d getScale(Vector3d dest) {
        return null;
    }

    @Override
    public Vector3d positiveZ(Vector3d dir) {
        return null;
    }

    @Override
    public Vector3d normalizedPositiveZ(Vector3d dir) {
        return null;
    }

    @Override
    public Vector3d positiveX(Vector3d dir) {
        return null;
    }

    @Override
    public Vector3d normalizedPositiveX(Vector3d dir) {
        return null;
    }

    @Override
    public Vector3d positiveY(Vector3d dir) {
        return null;
    }

    @Override
    public Vector3d normalizedPositiveY(Vector3d dir) {
        return null;
    }

    @Override
    public Matrix3d add(Matrix3dc other, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d sub(Matrix3dc subtrahend, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d mulComponentWise(Matrix3dc other, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d lerp(Matrix3dc other, double t, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateTowards(Vector3dc direction, Vector3dc up, Matrix3d dest) {
        return null;
    }

    @Override
    public Matrix3d rotateTowards(double dirX, double dirY, double dirZ, double upX, double upY, double upZ, Matrix3d dest) {
        return null;
    }

    @Override
    public Vector3d getEulerAnglesZYX(Vector3d dest) {
        return null;
    }

    @Override
    public Matrix3d obliqueZ(double a, double b, Matrix3d dest) {
        return null;
    }

    @Override
    public boolean equals(Matrix3dc m, double delta) {
        return false;
    }
}
