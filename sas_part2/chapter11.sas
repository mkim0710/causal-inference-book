/***************************************************************
PROGRAM 11.1
Sample averages by treatment level
Data from Figures 11.1 and 11.2
***************************************************************/;

data fig1;
input a y;
cards;
1 200
1 150
1 220
1 110
1 50
1 180
1 90
1 170
0 170
0 30
0 70
0 110
0 80
0 50
0 10
0 20
;
run;

proc plot data= fig1; 
	plot Y*A; 
run;

proc means data=fig1;
	class A;
	var Y;
run;


data fig2;
input a y;
cards;
1 110
1 80
1 50
1 40
2 170
2 30
2 70
2 50
3 110
3 50
3 180
3 130
4 200
4 150
4 220
4 210
;
run;

proc plot data= fig2; 
	plot Y*A; 
run;

proc means data=fig2;
	class A;
	var Y;
run;



/***************************************************************
PROGRAM 11.2
2-parameter linear model
Data from Figures 11.3 and 11.1
***************************************************************/;

data fig3;
input a y;
cards;
3   21	
11	54
17	33
23	101
29	85
37	65
41	157
53	120
67	111
79	200
83	140
97	220
60	230
71	217
15	11
45  190
;
proc plot data= fig3; 
	plot Y*A; 
run;

proc glm data= fig3;
	model Y= A/ clparm;
	estimate 'A=90' intercept 1 A 90;
run;

proc glm data= fig1;
	model Y= A/ clparm;
run;
quit;



/**********************************************www*****************
PROGRAM 11.3
3-parameter linear model
Data from Figure 11.3
***************************************************************/;

data fig4;
	set fig3;
	Asq= A*A;
run;

proc glm data= fig4;
	model Y = A Asq/ clparm;
	estimate 'A=90' intercept 1 A 90 Asq 8100;
run;
quit;
