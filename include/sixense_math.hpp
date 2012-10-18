/*
 *
 * SIXENSE CONFIDENTIAL
 *
 * Copyright (C) 2011 Sixense Entertainment Inc.
 * All Rights Reserved
 *
 */

#ifndef SIXENSE_MATH_HPP
#define SIXENSE_MATH_HPP

#include <math.h>
#include <iostream>
#include <string>

namespace sixenseMath {

	class Vector2 {
	public:
		Vector2();
		Vector2( const float x, const float y );
		float& operator [](const int idx);
		Vector2 operator -(const Vector2 rhs);
		Vector2 operator +(const Vector2 rhs);

		Vector2& operator +=(const Vector2& rhs);
		Vector2& operator *=(const float& rhs);
		Vector2& operator -=(const Vector2& rhs);
		float operator *(const Vector2 rhs); // dot product
		Vector2 operator *(const float rhs);
		Vector2 operator /(const float rhs);

		bool operator ==(const Vector2& rhs);

		void normalize();
		float length();
		void print( const std::string name=std::string() );
		void fill( float vec[2] );

	protected:
		float _vec[2];
	};

	class Vector3 {
	public:
		Vector3();
		Vector3( const Vector3& );
		Vector3( const float vec[3] );
		Vector3( const float x, const float y, const float z );
		float& operator [](const int idx);
		Vector3 operator -(const Vector3 rhs);
		Vector3 operator +(const Vector3 rhs);
		Vector3 operator ^(const Vector3 rhs); // cross product
		Vector3& operator +=(const Vector3& rhs);
		Vector3& operator *=(const float& rhs);
		Vector3& operator -=(const Vector3& rhs);
		float operator *(const Vector3 rhs); // dot product
		Vector3 operator *(const float rhs);
		Vector3 operator /(const float rhs);
		void normalize();
		float length();
		void print( const std::string name=std::string() );
		void fill( float vec[3] );

		bool operator ==(const Vector3& rhs);

		// Construction helpers
		static Vector3 normalize( const Vector3 );

	protected:
		float _vec[3];
	};

	class Vector4 {
	public:
		Vector4();
		Vector4( const Vector4& );
		Vector4( const Vector3&, float w );
		Vector4( const float vec[4] );
		Vector4( const float x, const float y, const float z, const float w );
		float& operator [](const int idx);
		Vector4 operator -(const Vector4 rhs) const;
		Vector4 operator ^(const Vector4 rhs) const; // cross product
		Vector4 operator +(const Vector4 rhs) const;
		float operator *(const Vector4 rhs); // dot product
		Vector4 operator *(const float rhs) const;
		Vector4 operator -(const float rhs) const;
		Vector4 operator /(const float rhs) const;

		bool operator ==(const Vector4& rhs);

		void normalize();
		float length();
		void print( const std::string name=std::string() );
		Vector4 operator *(const class Matrix4 rhs) const;
		void fill( float vec[4] );

	protected:
		float _vec[4];
	};

	class Quat : public Vector4 {
	public:
		Quat();
		Quat( const float x, const float y, const float z, const float w );
		Quat( const class Matrix3& );
		Quat( const Vector4& );
		Quat( const Vector3 xyz, float w );
		void print( const std::string name=std::string() );
		Vector3 operator *(const Vector3 rhs) const;
		Quat operator *(const Quat rhs) const;
		float dot(const Quat rhs) const; // dot product (operator * is already used for quat...)
		void invert(); // invert in place
		Quat inverse(); // leave this alone and return inverted copy
		Vector3 getEulerAngles();

		// Construction helpers
		static Quat rotation( const Vector3 from_vec, const Vector3 to_vec );
		static Quat rotation( const float angle_in_rad, const Vector3 axis );
		static Quat rotation( const Vector3 hpr_in_rad );
		static Quat rotation( const float heading, const float pitch, const float roll );
		static Quat slerp( const float t, const Quat a, const Quat b);
	};

	class Matrix3 {
	public:
		Matrix3();
		Matrix3( const Matrix3& );
		Matrix3( const float mat[3][3] );
		Matrix3( const float m00, const float m10, const float m20, const float m01, const float m11, const float m21, const float m02, const float m12, const float m22 );
		void fill( float mat[3][3] );
		Matrix3( const Vector3 col0, const Vector3 col1, const Vector3 col2 );
		Vector3& operator [](const int idx);
		Matrix3 operator *(const Matrix3 rhs);
		Matrix3 operator *(const float rhs);
		Matrix3 operator /(const float rhs);
		Matrix3 operator +(const Matrix3 rhs);
		Matrix3 operator -(const Matrix3 rhs);
		Matrix3 operator *(const Quat rhs);
		Vector3 operator *(const Vector3 rhs);
		Vector3 col( const int );
		Vector3 row( const int );
		void set_col( const int which, const Vector3 col );
		float trace();
		bool is_identity();
		void transpose();
		void print( const std::string name=std::string() );
		Vector3 getEulerAngles();

		// Construction helpers
		static Matrix3 rotation( const float angle_in_rad, const Vector3 axis );
		static Matrix3 rotation( const Vector3 hpr_in_rad );
		static Matrix3 rotation( const Quat rot );
		static Matrix3 rotation( const Vector3 from, const Vector3 to );
		static Matrix3 translation( const Vector3 trans );
		static Matrix3 scale( const float, const float, const float );
		static Matrix3 scale( const float );
		static Matrix3 transpose( const Matrix3 );

	protected:
		Vector3 _cols[3];
	};

	class Matrix4 {
	public:
		Matrix4();
		Matrix4( const Matrix4& );
		Matrix4( const Matrix3& );
		Matrix4( const float mat[4][4] );
		Matrix4( const float m00, const float m10, const float m20, const float m30, const float m01, const float m11, const float m21, const float m31, const float m02, const float m12, const float m22, const float m32, const float m03, const float m13, const float m23, const float m33 );
		void fill( float mat[4][4] );
		Matrix4( const Vector4 col0, const Vector4 col1, const Vector4 col2, const Vector4 col3 );
		Vector4& operator [](const int idx);
		Matrix4 operator *(const Matrix4 rhs);
		Matrix4 operator *(const float rhs);
		Matrix4 operator /(const float rhs);
		Matrix4 operator +(const Matrix4 rhs);
		Matrix4 operator -(const Matrix4 rhs);
		Matrix4 operator *(const Quat rhs);
		Vector4 operator *(const Vector4 rhs);
		Vector4 col( const int );
		Vector4 row( const int );
		void set_col( const int which, const Vector4 col );
		float trace();
		bool is_identity();
		void transpose();
		void print( const std::string name=std::string() );
		Vector3 getEulerAngles();

		// Construction helpers
		static Matrix4 rotation( const float angle_in_rad, const Vector3 axis );
		static Matrix4 rotation( const Quat rot );
		static Matrix4 rotation( const Vector3 from, const Vector3 to );
		static Matrix4 rotation( const Vector3 hpr_in_rad );
		static Matrix4 translation( const Vector3 trans );
		static Matrix4 scale( const float, const float, const float );
		static Matrix4 scale( const float );
		static Matrix4 transpose( const Matrix4 );

	protected:
		Vector4 _cols[4];
	};

	class Line {
		friend class Plane;
		public:
			Line( const Line& );
			Line( const Vector3& dir, const Vector3& pos ); 

			Vector3 getClosestPoint( const Vector3& );
			
		private:
			Vector3 _dir;
			Vector3 _pos1;
			Vector3 _pos2;
	};


	class Plane	{
	public:
		Plane(); 
		Plane( const Plane&	);
		Plane( Vector3 p0, Vector3 p1, Vector3 p2 );
		Plane( Vector3 point, Vector3 normal );

		void init();
		double whichSide( Vector3 p );
		Vector3 getClosestPoint( Vector3 in );
		Vector3 intersect( const Line line );
		Vector3 getNormal();

	private:
		double _a, _b, _c, _d;
		Vector3	_norm;
		Vector3	_p0, _p1, _p2;

	};
	
	inline float fix_angle(float angle,float center = 0) {
		float test_angle = angle;
		int cnt = 1;
		while ((test_angle-center) > M_PI) {
			test_angle = angle - cnt * 2*M_PI;
			cnt++;
		}
		angle = test_angle;
		cnt = 1;
		while ((test_angle-center) < -M_PI) {
			test_angle = angle + cnt * 2*M_PI;
			cnt++;
		}
		return test_angle;
	}

#ifndef TB_ANGLES_CLOSE_ENOUGH
#define TB_ANGLES_CLOSE_ENOUGH 0.0001
#endif

	class tb_angles {
	public:
		float yaw;
		float yaw_deg;
		float pitch;
		float pitch_deg;
		float roll;
		float roll_deg;
	
		tb_angles(float yaw, float pitch, float roll,bool deg=true) {
			if (deg) {
				this->yaw = yaw * M_PI / 180.;
				this->yaw_deg = yaw;
				this->pitch = pitch * M_PI / 180.;
				this->pitch_deg = pitch;
				this->roll = roll * M_PI / 180.;
				this->roll_deg = roll;
			} else {
				this->yaw = yaw;
				this->yaw_deg = yaw * 180. / M_PI;
				this->pitch = pitch;
				this->pitch_deg = pitch * 180. / M_PI;
				this->roll = roll;
				this->roll_deg = roll * 180. / M_PI;
			}
		}
		tb_angles(sixenseMath::Quat& q) { sixenseMath::Matrix3 m = sixenseMath::Matrix3::rotation(q); init(m); }
		tb_angles(sixenseMath::Matrix3& R) { init(R); }

		template<typename OtherType>
		tb_angles(OtherType R) {
			sixenseMath::Matrix3 ssR;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					ssR[i][j] = R[i][j];
				}
			}
			init(ssR);
		}
	
		sixenseMath::Quat toQuaternion() const { sixenseMath::Quat(toMatrix()); }
		sixenseMath::Matrix4 toTransform() const { return sixenseMath::Matrix4(toMatrix()); }
		sixenseMath::Matrix3 toMatrix() const {
			sixenseMath::Matrix3 Ryaw(
					cos(yaw), -sin(yaw), 0,
					sin(yaw),  cos(yaw), 0,
					0,         0,        1);
			sixenseMath::Matrix3 Rpitch(
					 cos(pitch), 0, sin(pitch),
					 0,          1, 0,
					-sin(pitch), 0, cos(pitch));
			sixenseMath::Matrix3 Rroll(
					1,  0,          0,
					0,  cos(roll), -sin(roll),
					0,  sin(roll),  cos(roll));
			return Ryaw * Rpitch * Rroll;
		}
	
	private:
		void init(sixenseMath::Matrix3& R) {
			yaw = 0;
			pitch = 0;
			roll = 0;

			bool skip = false;
			if (fabs(R[0][1]-R[1][0]) < TB_ANGLES_CLOSE_ENOUGH && fabs(R[0][2]-R[2][0]) < TB_ANGLES_CLOSE_ENOUGH && fabs(R[1][2]-R[2][1]) < TB_ANGLES_CLOSE_ENOUGH) {
				//matrix is symmetric
				if (fabs(R[0][1]+R[1][0]) < TB_ANGLES_CLOSE_ENOUGH && fabs(R[0][2]+R[2][0]) < TB_ANGLES_CLOSE_ENOUGH && fabs(R[1][2]+R[2][1]) < TB_ANGLES_CLOSE_ENOUGH) {
					//diagonal
					if (R[0][0] > 0) {
						if (R[1][1] > 0) {
							skip = true;
						} else {
							roll = M_PI;
						}
					} else if (R[1][1] > 0) {
						yaw = M_PI;
						pitch = M_PI;
					} else {
						yaw = M_PI;
					}
					skip = true;
				}
			}

			if (!skip) {
				sixenseMath::Vector3 vx = R * sixenseMath::Vector3(1,0,0);
				sixenseMath::Vector3 vy = R * sixenseMath::Vector3(0,1,0);

				yaw = atan2(vx[1],vx[0]);
				pitch = atan2(-vx[2], sqrt(vx[0]*vx[0] + vx[1]*vx[1]));

				sixenseMath::Matrix3 Ryaw(
							 cos(yaw), -sin(yaw), 0,
							 sin(yaw),  cos(yaw), 0,
							 0,         0,        1);
				sixenseMath::Matrix3 Rpitch(
						 cos(pitch), 0, sin(pitch),
						 0,          1, 0,
						-sin(pitch), 0, cos(pitch));
				sixenseMath::Vector3 vyp = Ryaw * Rpitch * sixenseMath::Vector3(0,1,0);
				sixenseMath::Vector3 vzp = Ryaw * Rpitch * sixenseMath::Vector3(0,0,1);

				float coeff = (vzp * vy) >= 0 ? 1 : -1;

				roll = coeff * acos(vyp *vy);
			}

			yaw_deg = yaw * 180. / M_PI;
			pitch_deg = pitch * 180. / M_PI;
			roll_deg = roll * 180. / M_PI;
		}
	};




	#include "sixense_math.cpp"

}


#endif
