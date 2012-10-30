Exercise 2

Given the vectors u = (x1*a,y1*b), v = (x2*a,y2*b), w = (x3*a,y3*b), 0 = (0,0)

1: <u , v> 			= <(x1*a,y1*b) , (x2*a,y2*b)>
					= x1*a*x2*a + y1*b*y2*b
					= x2*a*x1*a + y2*b*y1*b
					= <v,u>

2: <u + v , w>		= <(x1*a,y1*b)+(x2*a,y2*b) , (x3*a,y3*b)>
					= <(x1*a + x2*a , y1*b + y2*b) , (x3*a,y3*b)>
					= (x1*a + x2*a)*x3*a + (y1*a + y2*a)*y3*a
					= x1*a*x3*a + x2*a*x3*a + y1*a*y3*a + y2*a*y3*a
					= <u , w> + <v , w >

3: <c*u , v>		= <c*(x1*a , y1*b) , (x2*a,y2*b)>	
					= <(c*x1*a , c*y1*b) , (x2*a,y2*b)>	
					= c*x1*a*x2*a + c*y1*b*y2*b
					= c*( x1*a*x2*a + y1*b*y2*b ) 
					= c* <(x1*a , y1*b) , (x2*a,y2*b)>
					= c<u , v>

4: if x1*a > 0 and y1*b > 0

		<u , u > 	=  <(x1*a,y1*b) , (x1*a,y1*b)> 
					= x1*a*x1*a + y1*b*y2*b
					= (x1*a)**2 + (y1*b)**2 >= 0

   if x1*a > 0 and y1*b < 0

   		<u , u > 	=  <(x1*a,y1*b) , (x1*a,y1*b)> 
					= x1*a*x1*a + y1*b*y2*b
					= (x1*a)**2 + (-y1*b)**2			
					= (x1*a)**2 + (y1*b)**2 >= 0			

   if x1*a < 0 and y1*b < 0

   		<u , u > 	=  <(x1*a,y1*b) , (x1*a,y1*b)> 
					= x1*a*x1*a + y1*b*y2*b
					= (-x1*a)**2 + (-y1*b)**2			
					= (x1*a)**2 + (y1*b)**2 >= 0	

   if x1*a = 0 and y1*b = 0

   		<u , u > 	=  <(x1*a,y1*b) , (x1*a,y1*b)> 
					= x1*a*x1*a + y1*b*y2*b
					= (x1*a)**2 + (y1*b)**2
					= (0)**2 + (0)**2 = 0			
					