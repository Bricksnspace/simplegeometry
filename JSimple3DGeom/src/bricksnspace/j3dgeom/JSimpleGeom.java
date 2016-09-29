/*
	Copyright 2014 Mario Pascucci <mpascucci@gmail.com>
	This file is part of JSimple3DGeom

	JSimple3DGeom is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JSimple3DGeom is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with JSimple3DGeom.  If not, see <http://www.gnu.org/licenses/>.

*/


package bricksnspace.j3dgeom;



/**
 * Simple geometric computation 
 * 
 * @author Mario Pascucci
 *
 */
public class JSimpleGeom {

	
	/**
	 * Computes intersection between line and plane 
	 * Line is defined as two point in space
	 * Plane is defined as three point in space
	 * @param p plane defined by three points in a vector [x1,y1,z1,x2,y2,z2,...]
	 * @param l a line defined by two points in a vector [x1,y1,z1,x2,y2,z2]
	 * @return
	 */
	public static float[] planeXline(float[] p, Point3D l0, Point3D l1) {
		
		float a = l0.x - l1.x;	// xa-xb
		float b = p[3] - p[0];	// x1-x0
		float c = p[6] - p[0];	// x2-x0
		float d = l0.y - l1.y;	// ya-yb
		float e = p[4] - p[1];	// y1-y0
		float f = p[7] - p[1];	// y2-y1
		float g = l0.z - l1.z;	// za-zb
		float h = p[5] - p[2];	// z1-z0
		float i = p[8] - p[2];	// z2-z0
		float determinant = a*(e*i-f*h) - b*(i*d-f*g) + c*(d*h-e*g);
		if (determinant != 0.0f) {
			float A = (e*i-f*h)/determinant;
			float B = -(b*i-c*h)/determinant;
			float C = (b*f-c*e)/determinant;
			float t = A * (l0.x-p[0]) + B * (l0.y-p[1]) + C * (l0.z-p[2]);
			return new float[] {t,-a*t+l0.x,-d*t+l0.y,-g*t+l0.z};
		}
		return new float[] {-1f,0,0,0};
	}

	
	/**
	 * Computes intersection between line and plane 
	 * Line is defined as two point in space
	 * Plane is defined as three point in space
	 * @param p plane defined by three points in a vector [x1,y1,z1,x2,y2,z2,...]
	 * @param l a line defined by two points in a vector [x1,y1,z1,x2,y2,z2]
	 * @return
	 */
	public static float[] planeXline(Point3D p0, Point3D p1, Point3D p2, Point3D l0, Point3D l1) {
		
		float a = l0.x - l1.x;	// xa-xb
		float b = p1.x - p0.x;	// x1-x0
		float c = p2.x - p0.x;	// x2-x0
		float d = l0.y - l1.y;	// ya-yb
		float e = p1.y - p0.y;	// y1-y0
		float f = p2.y - p0.y;	// y2-y1
		float g = l0.z - l1.z;	// za-zb
		float h = p1.z - p0.z;	// z1-z0
		float i = p2.z - p0.z;	// z2-z0
		float determinant = a*(e*i-f*h) - b*(i*d-f*g) + c*(d*h-e*g);
		if (determinant != 0.0f) {
			float A = (e*i-f*h)/determinant;
			float B = -(b*i-c*h)/determinant;
			float C = (b*f-c*e)/determinant;
			float t = A * (l0.x-p0.x) + B * (l0.y-p0.y) + C * (l0.z-p0.z);
			return new float[] {t,-a*t+l0.x,-d*t+l0.y,-g*t+l0.z};
		}
		return new float[] {-1f,0,0,0};
	}

	
	/**
	 * returns distance^2 from point p to line a-b
	 * to computes minimum distance is not necessary to use square root
	 * you can use the square of distance 
	 * 
	 * @param p  point coordinates
	 * @param a first point of line
	 * @param b second point of line
	 * @return distance^2 from p to line a-b
	 */
	public static float point2line(Point3D p, Point3D a, Point3D b) {
		
		// cross product, vector ba by vector ap: (ba)x(ap) 
		// a = b-a  b = a-p
		float a1 = b.x-a.x;
		float a2 = b.y-a.y;
		float a3 = b.z-a.z;
		float b1 = a.x-p.x;
		float b2 = a.y-p.y;
		float b3 = a.z-p.z;
		float xp = a2*b3-a3*b2;
		float yp = a3*b1-a1*b3;
		float zp = a1*b2-a2*b1;
		float mod = a1*a1+a2*a2+a3*a3;
		if (mod == 0) {
			// degenerated line
			// return a big distance
			return 1.0e30f;
		}
		return (xp*xp+yp*yp+zp*zp) / mod;
	}
	
	
	/**
	 * return a normal vector (origin at 0,0,0) of two vector
	 * computed from three points (v1 = p2-p1, v2 = p3-p2)
	 * @param p1 first point
	 * @param p2 second point (common to two vectors)
	 * @param p3 third point
	 * @return a normal vector as Point3D (origin at 0,0,0) 
	 */
	public static Point3D calcNormal(Point3D p1, Point3D p2, Point3D p3) {

		float v1x,v1y,v1z,v2x,v2y,v2z,xn,yn,zn;

		v1x = p2.x - p1.x;
		v1y = p2.y - p1.y;
		v1z = p2.z - p1.z;
		v2x = p3.x - p2.x;
		v2y = p3.y - p2.y;
		v2z = p3.z - p2.z;
		xn = v1y * v2z - v1z * v2y;
		yn = v1z * v2x - v1x * v2z;
		zn = v1x * v2y - v1y * v2x;
		return new Point3D(xn,yn,zn);
	}
	
	
	/**
	 * returns cross product of vector p1-p2 by vector p3-p4 as 
	 * a vector originating on (0,0,0)
	 * 
	 * @param p1 first vector origin
	 * @param p2 first vector head
	 * @param p3 second vector origin
	 * @param p4 second vector head
	 * @return resulting vector head, originating at (0,0,0)  
	 */
	public static Point3D crossProd(Point3D p1, Point3D p2, Point3D p3, Point3D p4){

		float v1x,v1y,v1z,v2x,v2y,v2z,xn,yn,zn;

		v1x = p2.x - p1.x;
		v1y = p2.y - p1.y;
		v1z = p2.z - p1.z;
		v2x = p4.x - p3.x;
		v2y = p4.y - p3.y;
		v2z = p4.z - p3.z;
		xn = v1y * v2z - v1z * v2y;
		yn = v1z * v2x - v1x * v2z;
		zn = v1x * v2y - v1y * v2x;
		return new Point3D(xn,yn,zn);

	}
	
	
	/**
	 * return a cross product for two vector as a point3d
	 * @param v1
	 * @param v2
	 * @return
	 */
	public static Point3D crossVector(Point3D v1,Point3D v2) {
		
		return new Point3D(
				v1.y * v2.z - v1.z * v2.y,
				v1.z * v2.x - v1.x * v2.z,
				v1.x * v2.y - v1.y * v2.x);
	}
	
	
	/**
	 * return dot product of two vectors
	 * @param p1 first vector origin
	 * @param p2 first vector head
	 * @param p3 second vector origin
	 * @param p4 second vector head
	 * @return |v1|*|v2|*cos theta
	 */
	public static float dotPoints(Point3D p1, Point3D p2, Point3D p3, Point3D p4) {
		
		return (p2.x-p1.x)*(p4.x-p3.x)+(p2.y-p1.y)*(p4.y-p3.y)+(p2.z-p1.z)*(p4.z-p3.z);
	}
	
	
	public static float dotVector(Point3D p1, Point3D p2) {
		
		return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
	}
	
	
	
	/**
	 * modulo of a vector
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static float modulo(Point3D p1, Point3D p2) {
		
		float x = p2.x-p1.x; 
		float y = p2.y-p1.y;
		float z = p2.z-p1.z;
		return (float)Math.sqrt(x*x+y*y+z*z);
	}
	
	
	/**
	 * return angle between vectors p2-p1/p4-p3 as a cos(theta) 
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param p4
	 * @return
	 */
	public static float dotProdAngle(Point3D p1, Point3D p2, Point3D p3, Point3D p4) {
		
		float mod = modulo(p1, p2)*modulo(p3, p4);
		if (mod == 0) 
			return 1;
		return dotPoints(p1, p2, p3, p4)/(mod);
	}
	
	/**
	 * Distance between two segments, from:
	 * http://mathworld.wolfram.com/Line-LineDistance.html
	 * 
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param p4
	 * @return
	 */
	public static float line2line(Point3D p1, Point3D p2, Point3D p3, Point3D p4) {
		
		Point3D cross = crossProd(p1, p2, p3, p4);
		if (cross.getDistSq(new Point3D(0, 0, 0)) == 0) {
			// parallel lines
			return (float) Math.sqrt(p1.getDistSq(p3));
		}
		Point3D c = p1.vector(p3);
		return (float) Math.abs(dotVector(c, cross))/cross.modulo();
	}
	
	
	/**
	 * return a point on segment s1-s2 that is nearest to line l1-l2
	 * from: http://stackoverflow.com/a/18994296/2772982
	 * @return
	 */
	public static Point3D nearestPointToSegments(Point3D s1, Point3D s2, Point3D l1, Point3D l2, boolean clamp) {
		
		Point3D a = s1.vector(s2);
		Point3D b = l1.vector(l2);
		Point3D _a = a.normalize();
		Point3D _b = b.normalize();
		
		Point3D cross = crossVector(_a, _b);
		float denom = cross.getDistSq(new Point3D(0,0,0));
		// parallel lines? return s1 as a nearest point
		if (denom == 0) {
			return s1;
		}
		Point3D t = s1.vector(l1);
		// matrix 3x3 given by 
		//   tx  ty  tz
		//  _bx _by _bz
		//  crx cry crz
		//det == a*(e*i-f*h) - b*(i*d-f*g) + c*(d*h-e*g);
		float det = t.x*(_b.y*cross.z-_b.z*cross.y) - 
				t.y*(_b.x*cross.z-_b.z*cross.x) +
				t.z*(_b.x*cross.y-_b.y*cross.x);
		float t0 = det/denom;
		if (clamp) {
			if (t0 < 0)
				return s1;
			if (t0 > a.modulo())
				return s2;
		}
		return new Point3D(s1.x+(_a.x*t0),s1.y+(_a.y*t0),s1.z+(_a.z*t0));
	}

	
	
	/**
	 * Computes a rotation matrix around vector cp1-cp2 by angle degree
	 * @param cp1 base point of vector
	 * @param cp2 head point of vector
	 * @param angle rotation angle in degree
	 * @return rotation matrix
	 */
	public static Matrix3D axisRotMatrix(Point3D cp1, Point3D cp2, float angle) {
	
		// from http://www.gamedev.net/page/resources/_/technical/math-and-physics/do-we-really-need-quaternions-r1199
		// http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
		Point3D r = cp1.vector(cp2).normalize();
		float rad = (float) (angle*Math.PI/180);
		float c = (float) Math.cos(rad);
		float s = (float) Math.sin(rad);
		float t = 1 - c;
		return new Matrix3D(
				t*r.x*r.x+c,	t*r.x*r.y-s*r.z,	t*r.x*r.z+s*r.y,
				t*r.x*r.y+s*r.z,t*r.y*r.y+c,		t*r.y*r.z-s*r.x,
				t*r.x*r.z-s*r.y,t*r.y*r.z+s*r.x,	t*r.z*r.z+c,
				0,0,0
				);
	}



	

	/**
	 * return a transformation matrix to align an object to a plane defined from three points 
	 * 
	 * @param p1 first point 
	 * @param p2 second point
	 * @param p3 third point
	 * @return a transformation matrix needed to align an object to given plane
	 */
	public static Matrix3D alignMatrix(Point3D p1, Point3D p2, Point3D p3) {
		
		float angley, anglez;
		Point3D normal = calcNormal(p1, p2, p3);
		if (normal.z == 0 && normal.x == 0) {
			angley = 0;
			if (Math.signum(normal.y) < 0f) {
				return Matrix3D.getRotateZ((float)-Math.PI).moveTo(p1.x, p1.y, p1.z);
			}
			else {
				return new Matrix3D().moveTo(p1.x, p1.y, p1.z);
			}
		}
		// compute rotations
		angley = (float)Math.atan2(normal.z, normal.x);
		Matrix3D m = Matrix3D.getRotateY(-angley);
		Point3D projxy = new Point3D(m.transformNormal(normal.x, normal.y, normal.z));
		anglez = (float)Math.atan2(projxy.x,projxy.y);
		return Matrix3D.getRotateZ(anglez).rotateY(angley).moveTo(p1.x, p1.y, p1.z);
	}


//	/**
//	 * return a transformation matrix to align an object aligned to Z axis 
//	 * to a vector defined by point p1-p2 
//	 * @param p1 base of vector
//	 * @param p2 head of vector
//	 * @return a transformation matrix needed to align to given vector
//	 */
//	public static Matrix3D alignZvectorMatrix(Point3D p1, Point3D p2) {
//		
//		float angley, anglex;
//		Point3D normal = p1.vector(p2).normalize();
//		if (normal.z <= 0.01 && normal.x <= 0.01) {
//			float anglez = (float) Math.atan(normal.x/normal.y);
//			Matrix3D m = Matrix3D.getRotateZ(-anglez);
//			Point3D projyz = new Point3D(m.transformNormal(normal.x, normal.y, normal.z));
//			anglex = (float)Math.atan2(projyz.z,projyz.y);
//			return Matrix3D.rotX90().rotateX(-anglex).rotateZ(anglez);
//		}
//		// compute rotations
//		if (normal.z == 0)
//			angley = (float) Math.PI/2;
//		else 
//			angley = (float)Math.atan(normal.x/normal.z);
//		Matrix3D m = Matrix3D.getRotateY(angley);
//		Point3D projxy = new Point3D(m.transformNormal(normal.x, normal.y, normal.z));
//		anglex = (float)Math.atan2(projxy.y,projxy.z);
//		//System.out.println("z="+anglez+" y="+angley);
//		return Matrix3D.getRotateX(anglex).rotateY(-angley);
//	}


	/** 
	 * rotation matrix to align two vectors in same direction and orientation
	 * @param cpx vector to align
	 * @param tpx target vector
	 * @return
	 */
	public static Matrix3D alignMatrix(Point3D cp1, Point3D cp2, Point3D tp1, Point3D tp2) {
	
		// from http://www.gamedev.net/page/resources/_/technical/math-and-physics/do-we-really-need-quaternions-r1199
		// http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
		Point3D vx = cp1.vector(cp2).normalize();//cp.getP1().vector(cp.getP2()).normalize();
		Point3D vy = tp1.vector(tp2).normalize();//to.getP1().vector(to.getP2()).normalize();
		Point3D r = crossVector(vx, vy).normalize();
		float c = dotVector(vx, vy);
		// rounding error and precision loss -> c == 1.00000001
		if (c > 1)
			c = 1;
		float s = (float) Math.sin(Math.acos(c));
		float t = 1 - c;
		return new Matrix3D(
				t*r.x*r.x+c,	t*r.x*r.y-s*r.z,	t*r.x*r.z+s*r.y,
				t*r.x*r.y+s*r.z,t*r.y*r.y+c,		t*r.y*r.z-s*r.x,
				t*r.x*r.z-s*r.y,t*r.y*r.z+s*r.x,	t*r.z*r.z+c,
				0,0,0
				);
	}


	/**
	 * minimum rotation to align (parallel) two vectors, ignoring orientation and 
	 * inverting first vector if needed
	 * @param cp1 vector origin to align
	 * @param cp2 vector direction
	 * @param tp1 target vector origin
	 * @param tp2 target vector direction
	 * @return
	 */
	public static Matrix3D minAlignMatrix(Point3D cp1, Point3D cp2, Point3D tp1, Point3D tp2) {
		
		// from http://www.gamedev.net/page/resources/_/technical/math-and-physics/do-we-really-need-quaternions-r1199
		// http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
		Point3D vx = cp1.vector(cp2).normalize();//cp.getP1().vector(cp.getP2()).normalize();
		Point3D vy = tp1.vector(tp2).normalize();//to.getP1().vector(to.getP2()).normalize();
		float c = dotVector(vx, vy);
		// invert a vector to get a littler angle
		if (c < 0) {
			vx = cp2.vector(cp1).normalize();//cp.getP2().vector(cp.getP1()).normalize();
			c = dotVector(vx, vy);
		}
		Point3D r = crossVector(vx, vy).normalize();
		float s = (float) Math.sin(Math.acos(c));
		float t = 1 - c;
		return new Matrix3D(
				t*r.x*r.x+c,	t*r.x*r.y-s*r.z,	t*r.x*r.z+s*r.y,
				t*r.x*r.y+s*r.z,t*r.y*r.y+c,		t*r.y*r.z-s*r.x,
				t*r.x*r.z-s*r.y,t*r.y*r.z+s*r.x,	t*r.z*r.z+c,
				0,0,0
				);
	}
	
	
	
	/**
	 * Returns a point on a cubic Bézier curve with four control points
	 * <p>
	 * from: http://it.wikipedia.org/wiki/Curva_di_B%C3%A9zier#Applicazioni_nella_computer_grafica
	 * @param cp0 initial control point 
	 * @param cp1 first control point
	 * @param cp2 second control point
	 * @param cp3 final control point
	 * @param t point on curve (0<=t<=1)
	 * @return point coordinates
	 */
	public static float[] pointOnCubicBezier(Point3D cp0, Point3D cp1, Point3D cp2, Point3D cp3, float t) {
		float   ax, bx, cx;
	 	float   ay, by, cy;
	 	float   az, bz, cz;
	 	float   tSquared, tCubed;
	 
	 	cx = 3.0f * (cp1.x - cp0.x);
	 	bx = 3.0f * (cp2.x - cp1.x) - cx;
	 	ax = cp3.x - cp0.x - cx - bx;
	 
	 	cy = 3.0f * (cp1.y - cp0.y);
	 	by = 3.0f * (cp2.y - cp1.y) - cy;
	 	ay = cp3.y - cp0.y - cy - by;
	 
	 	cz = 3.0f * (cp1.z - cp0.z);
	 	bz = 3.0f * (cp2.z - cp1.z) - cz;
	 	az = cp3.z - cp0.z - cz - bz;
	 
	 	// compute p(t)
	 	tSquared = t * t;
	 	tCubed = tSquared * t;
	 
	 	return new float[] {(ax * tCubed) + (bx * tSquared) + (cx * t) + cp0.x,
	 				(ay * tCubed) + (by * tSquared) + (cy * t) + cp0.y,
	 				(az * tCubed) + (bz * tSquared) + (cz * t) + cp0.z};
	}
	 
	
	
	/**
	 * generates a cubic Bézier curve from four control points
	 * @param cp0 first control point (start)
	 * @param cp1 second control point
	 * @param cp2 third control point
	 * @param cp3 fourth control point (end)
	 * @param numpoint how many point (start and end included)
	 * @param buffer user-supplied buffer of floats to hold generated points
	 * @param offset offset in point in buffer
	 */
	public static void generateBezier(Point3D cp0, Point3D cp1, Point3D cp2, Point3D cp3, 
			int numpoint, float[] buffer, int offset) {
		
		float   ax, bx, cx;
	 	float   ay, by, cy;
	 	float   az, bz, cz;
	 	float   tSquared, tCubed;
	 
		float step = 1f/(numpoint-1);
		float t = 0;
	 	cx = 3.0f * (cp1.x - cp0.x);
	 	bx = 3.0f * (cp2.x - cp1.x) - cx;
	 	ax = cp3.x - cp0.x - cx - bx;
	 
	 	cy = 3.0f * (cp1.y - cp0.y);
	 	by = 3.0f * (cp2.y - cp1.y) - cy;
	 	ay = cp3.y - cp0.y - cy - by;
	 
	 	cz = 3.0f * (cp1.z - cp0.z);
	 	bz = 3.0f * (cp2.z - cp1.z) - cz;
	 	az = cp3.z - cp0.z - cz - bz;
	 	offset *= 3;
		for (int i=0;i<numpoint*3;i+=3) {
		 	tSquared = t * t;
		 	tCubed = tSquared * t;
			buffer[offset+i] = (ax * tCubed) + (bx * tSquared) + (cx * t) + cp0.x;
			buffer[offset+i+1] = (ay * tCubed) + (by * tSquared) + (cy * t) + cp0.y;
			buffer[offset+i+2] = (az * tCubed) + (bz * tSquared) + (cz * t) + cp0.z;
			t += step;
		}
		buffer[offset+numpoint*3] = cp3.x;
		buffer[offset+numpoint*3+1] = cp3.y;
		buffer[offset+numpoint*3+2] = cp3.z;
	}
	
	
	
}
