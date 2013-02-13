#ifndef AUXILIARY_H
#define AUXILIARY_H

typedef double real;

enum Axis {XX, YY, ZZ};
enum Direction {BACKWARD, FORWARD};
enum Side {LEFT, BACK, BOTTOM, RIGHT, FRONT, TOP};

inline Side side (Axis axis, Direction direction)
{
	if (axis == XX) return direction == BACKWARD ? LEFT : RIGHT; 
	if (axis == YY) return direction == BACKWARD ? FRONT : BACK;
	else return direction == BACKWARD ? BOTTOM : TOP;
}
inline Side side (int a, int d) { return side (Axis (a), Direction (d)); }

inline Axis axis (Side side) { return Axis (side%3); }
inline Axis axis (int s) { return Axis (s%3); }

inline Direction direction (Side side) { return Direction (side%2); }
inline Direction direction (int s) { return Direction (s%2); }

inline Side reflect (Side s) { return side (axis (s), !direction (s)); }

inline real sqr(real x) { return x*x; }

#endif
