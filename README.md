SymKalkulator
=============

SymKalkulator is a desktop calculator app. It uses RPN (Reverse Polish Notation), fractions, symbols like π and φ, and physical units.
It is called "Kalkulator" rather than "calculator" because it is in Norwegian.

## Commands
The following commands are implemented:
```
Rydd stabelen                                  →  rydd 
Slett y fra stabelen                           → slett 
Bytt x og y                                    →  bytt 
Kopier y                                       →  kopi 
Adder x + y                                    →   +   
Subtraher x - y                                →   -   
Bytt fortegn av y                              →   --  
Multipliser x·y                                →   *   
Divider x/y                                    →   /   
Beregn y prosent av x                          →   %   
Beregn fakultet av y                           →   !   
Beregn binominalkoeffisient                    →   ()  
Beregn potens xʸ                               →   **  
Beregn kvadratrot av y                         →   v   
Beregn y-te rot av x                           →  rot  
Beregn briggsk logaritme av y                  →   lg  
Beregn naturlig logaritme av y                 →   ln  
Beregn binær logaritme y                       →   lb  
Beregn sinus av y                              →  sin  
Beregn cosinus av y                            →  cos  
Beregn tangens av y                            →  tan  
Beregn arcsinus av y                           → arcsin
Beregn arccosinus av y                         → arccos
Beregn arctangens av y                         → arctan
Beregn det y-te fibonacci-tallet               →  fib  
Beregn 1/y                                     →   //  
Beregn 1/y                                     → resiprok
Beregn rest av divisjonen x/y                  →  mod  
Beregn rest av divisjonen x/y                  →  rest 
Beregn summen av alle tallene i stabelen       →  sum  
Beregn gjennomsnitt av alle tallene i stabelen →   Ø   
Beregn produktet av alle tallene i stabelen    →  prod 
Omgjør y til heltall                           →  hel  
Omgjør y til enheten <enhet>                   → <enhet>
Omgjør y til SI-enheter                        →   SI  
Utvid et uttrykk                               → utvid 
Evaluer et uttrykk                             →  eval 
Eulers tall e                                  →   e   
sirkeltallet π                                 →   π   
sirkeltallet π                                 →   pi  
Det dobbelte sirkeltallet τ                    →  tau  
Det irrasjonale tallet φ                       →   fi  
Det irrasjonale tallet φ                       →  phi  
uendelig ∞                                     →   oo  
gravitasjonskonstanten G                       →   G   
lysets hastighet c                             →   c   
Plancks redusete konstant ħ                    →  hbar 
Plancks redusete konstant ħ                    →   ħ   
Coulombs konstant kₑ                           →   k   
gasskonstanten R                               →   R   

```
If you enter a physical unit, the app checks whether the current number already has a unit. If the number has no unit, the unit is applied to the number. If the number has a unit, it is converted to the new unit. A warning occurs if the units are incompatible. Numbers that are multiples of π are treated as radians.

To enter a number together with a unit, you must type two spaces between the number and the unit. The app understands all units known to Sympy, and additionally `min` and `'` for minute, `"` for second, `a` for year, `d` for day, `grad` for degree, `in` for inch and `Å` for Ångstrøm.

## Example
The following screenshots show how you can calculate the volume of a cone with a base radius of 2.5 and a height of 3.8:

Type in `1/3`:![Screenshot of SymKalkulator. 1/3 is typed in.](/screenshots/Skjermbilde%202025-10-02%20113848.png)
Type in `pi`:![Screenshot of SymKalkulator. Pi is is typed in.](/screenshots/Skjermbilde%202025-10-02%20114033.png)
Type in `*`:![Screenshot of SymKalkulator. * is typed in to multiply 1/3 and pi.](/screenshots/Skjermbilde%202025-10-02%20114056.png)
Type in `2,5`:![Screenshot of SymKalkulator. 2,5 is typed in.](/screenshots/Skjermbilde%202025-10-02%20114143.png)
Type in `2**`:![Screenshot of SymKalkulator. 2** is typed in to square 2.5.](/screenshots/Skjermbilde%202025-10-02%20114203.png)
Type in `*`:![Screenshot of SymKalkulator. * is typed in to multiply pi/3 and 6.25.](/screenshots/Skjermbilde%202025-10-02%20114253.png)
Type in `3,8*`:![Screenshot of SymKalkulator. 3,8* is typed in to multiply 25pi/12 with 3.8.](/screenshots/Skjermbilde%202025-10-02%20114313.png)
The result is exactly 95π/12 or 24.87094184091919647116259345096273... with an error of 1e-32:![Screenshot of SymKalkulator. The result is 95pi/12 or 24.87094184091919647116259345096273... with an error of 1e-32.](/screenshots/Skjermbilde%202025-10-02%20114326.png)

**This cross-platform app was generated by [Briefcase](https://briefcase.readthedocs.io/) - part of [The BeeWare Project](https://beeware.org/).**
