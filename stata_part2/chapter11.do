/***************************************************************
PROGRAM 11.1 for Stata 11
Sample averages by treatment level
Data from Figures 11.1 and 11.2
***************************************************************/
clear
input A Y
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
end

scatter Y A, ylab(0(50)250) xlab(0 1) xscale(range(-0.5 1.5))
bysort A: sum Y
save fig1, replace
clear

input A Y
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
end

scatter Y A, ylab(0(50)250) xlab(0(1)4) xscale(range(0 4.5))
bysort A: sum Y
save fig2, replace
clear

/***************************************************************
PROGRAM 11.2 for Stata 11
2-parameter linear model
Data from Figures 11.3 and 11.1
***************************************************************/

input A Y
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
end

scatter Y A, ylab(0(50)250) xlab(0(10)100) xscale(range(0 100))
regress Y A, noheader cformat(%5.2f)
lincom _b[_cons] + 90*_b[A]

scatter Y A, ylab(0(50)250) xlab(0(10)100) xscale(range(0 100)) || lfit Y A

save fig3, replace
clear

use fig1, clear
regress Y A, noheader cformat(%5.2f)

/**********************************************www*****************
PROGRAM 11.3 for Stata 11
3-parameter linear model
Data from Figure 11.3
***************************************************************/

use fig3, clear
gen Asq = A*A
regress Y A Asq, noheader cformat(%5.2f)
lincom _b[_cons] + 90*_b[A] + 90*90*_b[Asq]

scatter Y A, ylab(0(50)250) xlab(0(10)100) xscale(range(0 100)) || qfit Y A
clear

