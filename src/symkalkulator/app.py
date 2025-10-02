"""
En kalkulator-app som bruker omvendt polsk notasjon (OPN) og symboler
"""

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER

# Matplotlib, PIL.Image og io brukes til å vise resultatet med LaTeX
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import PIL.Image
import io

# SymPy brukes for symbolisk matematikk
from sympy import *
from sympy.ntheory.digits import digits
# from sympy.core.numbers import equal_valued

import re # Regular Expressions for å analysere det som blir tastet inn
import pyperclip # for å kopiere verdier til utklippstavla

MAKS_SIFRE = Integer(32) # Maksimalt antall sifre som kan vises i resultatfeltet

# Stringsymboler
IKKETALL = u'\u22A5'
UENDELIG = u'\u221E'
PI = u'\u03C0'
TAU = u'\u03C4'
FI = u'\u03C6'

# Lage hjelpetekst
kommandoer = ['slett', 'bytt', 'kopi', '+', '-', '--', '*', '/', '%', '!', '()', '**',  'v', 'rot',
              'lg', 'ln', 'lb', 'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'fib', '//', 
              'mod', 'hel', 'grad', 'rad', 'utvid', 'eval', 'pi', 'tau', 'e', 'fi']
forklaringer = ['Slett y fra stabelen',
                'Bytt x og y',
                'Kopier y',
                'Adder x + y',
                'Subtraher x - y',
                'Bytt fortegn av y',
                'Multipliser x·y',
                'Divider x/y',
                'Beregn y prosent av x',
                'Beregn fakultet av y',
                'Beregn binominalkoeffisient',
                'Potenser xʸ',
                'Beregn kvadratrot av y',
                'Beregn y-te rot av x',
                'Beregn briggsk logaritme av y',
                'Beregn naturlig logaritme av y',
                'Beregn binær logaritme y',
                'Beregn sinus av y',
                'Beregn cosinus av y',
                'Beregn tangens av y',
                'Beregn arcsinus av y',
                'Beregn arccosinus av y',
                'Beregn arctangens av y',
                'Beregn det y-te fibonacci-tallet',
                'Beregn 1/y',
                'Beregn rest av divisjonen x/y',
                'Omgjør y til heltall',
                'Omgjør y fra radianer til grader',
                'Omgjør y fra grader til radianer',
                'Utvid et uttrykk',
                'Evaluer et uttrykk',
                'Det irrasjonale tallet ' + PI,
                'Det irrasjonale tallet ' + TAU,
                'Det irrasjonale tallet e',
                'Det irrasjonale tallet ' + FI]

hjelpetittel = 'Hvis x og y er nederste tallene, for <funksjon> tast <innput>\n\n'
hjelpetekst = ''
for kommando, funksjon in zip(kommandoer, forklaringer):
    hjelpetekst += f'{funksjon:<33} \u2192 {kommando:^6}\n'

varseltekst = '' # Global variable for varsler


def sjekk_resultat(resultat: Expr) -> tuple[Rational, Expr]:
    """Denne funksjonen sjekker resultatet av en beregning og returnerer det rasjonale og irrasjonale delen av resultatet.
       Hvis resultatet er ugyldig, returneres det som NaN og en varseltekst blir satt.
    Args:
        resultat (Expr): Tallet som skal sjekkes.

    Returns:
        Rational:   Rasjonalt del av resultatet.
        Expr:       Irrasjonalt del av resultatet.
    """
    global varseltekst

    forenklet_resultat = resultat.simplify()
    rasjonal, irrasjonal = forenklet_resultat.as_content_primitive()

    if irrasjonal == S.NaN:
        varseltekst = 'OBS! Ugyldig input. Resultatet er ikke et reelt tall.'
        rasjonal = S.NaN
        irrasjonal = S.One
        return rasjonal, irrasjonal
    if irrasjonal == S.Zero:
        rasjonal = S.Zero
        irrasjonal = S.One
        return rasjonal, irrasjonal
    if irrasjonal.is_extended_negative:
        rasjonal *= -1
        irrasjonal *= -1
    if irrasjonal == S.Infinity:
        rasjonal *= S.Infinity
        irrasjonal = S.One
    return rasjonal, irrasjonal


class tall:
    """Denne klassen er datastrukturen som brukes til å representere tall i SymKalkulator.
       Jeg har valgt å ikke bruke en utvidelse som representerer brøktall, men å lage min egen.
       Slik får jeg de funksjonene slik som jeg ønsker dem.
    """    
    def __init__(self, tall: str) -> None:
        """Denne funksjonen skaper et nytt tall-objekt ved å ta en string og sjekke hvilket tall den representerer.
           Hvert tall representeres som brøk med teller og nevner og en bool som angir om tallet er negativt.

        Args:
            tall (str): Tegnene som skal bli et tall, f.eks. '1', '0,5', '1/2', '1 2/3' eller '1,5e-3'
        """    
        self.rasjonal: Rational = S.NaN   # Den rasjonale delen av tallet
        self.irrasjonal: Expr = S.One  # Den irrasjonale delen av tallet (Det er 1 hvis tallet er rasjonalt)
        self.feil: Integer = S.Zero  # Feil i beregningen
        negativ: bool = False  # Er tallet negativt?

        if tall == '': # Det er ikke et tall
            return
        
        if tall[0] == '-':
            negativ = True
            tall = tall[1:]
        elif tall[0] == '+':
            tall = tall[1:]
        
        if tall == '' or tall == '0/0': # Det er ikke et tall
            return
        
        # Det fungerer dessverre ikke med match case. Derfor må jeg bruke if.
        if tall == '0':
            self.rasjonal = S.Zero
            return
        if tall == '1':
            self.rasjonal = S.NegativeOne if negativ else S.One
            return
        if tall == '0,5' or tall == '1/2':
            self.rasjonal = S.Half
            if negativ:
                self.rasjonal *= S.NegativeOne
            return
        if tall == UENDELIG or tall == '1/0' or tall == 'oo':
            self.rasjonal = S.NegativeInfinity if negativ else S.Infinity
            return
        if tall == PI or tall == 'pi':
            self.rasjonal = S.NegativeOne if negativ else S.One
            self.irrasjonal = S.Pi
            return
        if tall == TAU or tall == 'tau': # τ = 2π
            self.rasjonal = Integer(-2) if negativ else Integer(2)
            self.irrasjonal = S.Pi
            return
        if tall == E or tall == 'e':
            self.rasjonal = S.NegativeOne if negativ else S.One
            self.irrasjonal = S.Exp1
            return
        if tall == FI or tall == 'fi' or tall == 'phi': # Det gyldne snitt φ
            self.rasjonal = S.NegativeOne if negativ else S.One
            self.irrasjonal = S.GoldenRatio
            return
        if re.fullmatch(r'^\d+$', tall): # heltall
            self.rasjonal = Integer(tall)
            if negativ:
                self.rasjonal *= -1
            return
        if re.fullmatch(r'^\d+ \d+/\d+$', tall): # blandet tall
            hel, brøkdel = tall.split(' ')
            self.rasjonal = Rational(Integer(hel) + Rational(brøkdel))
            if negativ:
                self.rasjonal *= S.NegativeOne
            return
        if re.fullmatch(r'^\d+/\d+$', tall): # brøk
            teller, nevner = tall.split('/')
            self.rasjonal = Rational(teller, nevner)
            if negativ:
                self.rasjonal *= S.NegativeOne
            return
        if re.fullmatch(r'^\d+,\d+$', tall): # desimaltall
            self.rasjonal = Rational(tall.replace(',', '.'))
            if negativ:
                self.rasjonal *= S.NegativeOne
            return
        if re.fullmatch(r'^\d+(,\d+)?e[+-]?\d+$', tall): # vitenskapelig format
            self.rasjonal = Rational(tall.replace(',', '.'))
            if negativ:
                self.rasjonal *= S.NegativeOne
            return


    def kopi(self):
        """Lager en kopi av et tall

        Returns:
            tall: et nytt tall med samme verdi som self
        """        
        c = tall('1')
        c.rasjonal = self.rasjonal
        c.irrasjonal = self.irrasjonal
        return c


    def gjør_hel(self) -> None:
        """Tar bare heltallig delen av et tall, dvs. 2,5 blir 2 og -2,5 blir -2
        """
        if self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity: # ∞ eller -∞
            return
        if self.rasjonal != S.NaN:
            self.rasjonal = Integer(Mul(self.rasjonal, self.irrasjonal))
            self.irrasjonal = S.One 


    def resiprok(self) -> None:
        """Beregner resiprokverdien av et tall, dvs. 1/self
        """        
        if self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity: # 1/∞ = 0
            self.rasjonal = S.Zero
            return
        if self.rasjonal == S.Zero: # 1/0 = ∞ 
            self.rasjonal = S.Infinity
            return
        if self.rasjonal == S.One and self.irrasjonal == S.One: # 1/1 = 1
            return
        if self.rasjonal == S.NegativeOne and self.irrasjonal == S.One: # 1/1 = 1
            return
        if self.rasjonal != S.NaN:
            b = self.kopi()
            self.rasjonal = S.One/b.rasjonal
            self.irrasjonal = S.One/b.irrasjonal
            return


    def pluss(self, addend) -> None:
        """Plusser sammen to tall: self + addend

        Args:
            addend (tall): Tallet som legges til
        """
        global varseltekst
        if addend.__class__ != tall: # Kan bare addere et tall
            return
        if addend.rasjonal == S.NaN: # kan bare addere tall
            return
        if self.rasjonal == S.NaN: # NaN + x = x
            self.rasjonal = addend.rasjonal
            self.irrasjonal = addend.irrasjonal
            return
        if self.rasjonal == S.NegativeInfinity and addend.rasjonal == S.Infinity: # -∞ + ∞ er ikke definert
            varseltekst = '-' + UENDELIG + ' pluss ' + UENDELIG + ' er ikke definert.'
            self.rasjonal = S.NaN
            return
        if self.rasjonal == S.Infinity and addend.rasjonal == S.NegativeInfinity: # ∞ - ∞ er ikke definert
            varseltekst = UENDELIG + ' minus ' + UENDELIG + ' er ikke definert.'
            self.rasjonal = S.NaN
            return
        if self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity: # ∞ + x = ∞
            return
        if addend.rasjonal == S.Infinity or addend.rasjonal == S.NegativeInfinity: # x + ∞ = ∞
            self.rasjonal = addend.rasjonal
            self.irrasjonal = addend.irrasjonal
            return
        if self.irrasjonal == S.One and addend.irrasjonal == S.One: # Det er to rasjonale tall
            self.rasjonal += addend.rasjonal
            return

        verdi = Add(Mul(self.rasjonal, self.irrasjonal), Mul(addend.rasjonal, addend.irrasjonal))
        self.rasjonal, self.irrasjonal = sjekk_resultat(verdi)

    
    def minus(self, minuend) -> None:
        """Trekker fra et tall: self - minuend

        Args:
            minuend (tall): Tallet som trekkes fra
        """
        if minuend.__class__ != tall: # Kan bare subtrahere et tall
            return
        b = minuend.kopi()
        b.rasjonal *= -1 # Bytter fortegn på minuend
        self.pluss(b)

            
    def ganger(self, faktor) -> None:
        """Ganger to tall med hverandre: self * faktor

        Args:
            faktor (tall): Tallet som ganges med
        """        
        global varseltekst
        if faktor.__class__ != tall: # Kan bare gange tall
            return
        if self.rasjonal == S.NaN: # kan bare gange tall
            varseltekst = 'Du kan bare gange tall.'
            return
        if faktor.rasjonal == S.NaN: # kan bare gange med tall
            varseltekst = 'Du kan bare gange tall.'
            return
        if (self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity) \
            and (faktor.rasjonal == S.Infinity or faktor.rasjonal == S.NegativeInfinity): # ∞ * ∞ er ikke definert
            varseltekst = UENDELIG + ' ganger ' + UENDELIG + ' er ikke definert.'
            self.rasjonal = S.NaN
            return
        if self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity: # ∞ * x er ∞
            match faktor.rasjonal:
                case S.Zero:
                    varseltekst = UENDELIG + ' ganger 0 er ikke definert.'
                    self.rasjonal = S.NaN
                case S.Infinity:
                    varseltekst = UENDELIG + ' ganger ' + UENDELIG + ' er ikke definert.'
                    self.rasjonal = S.NaN
                case S.NegativeInfinity:
                    varseltekst = UENDELIG + ' ganger ' + UENDELIG + ' er ikke definert.'
                    self.rasjonal = S.NaN
                case _:
                    self.rasjonal = faktor.rasjonal * self.rasjonal
            self.irrasjonal = S.One
            return
        if faktor.rasjonal == S.Infinity or faktor.rasjonal == S.NegativeInfinity: # x * ∞ er ∞
            match self.rasjonal:
                case S.Zero:
                    varseltekst = '0 ganger ' + UENDELIG + ' er ikke definert.'
                    self.rasjonal = S.NaN
                case _:
                    self.rasjonal = faktor.rasjonal * self.rasjonal
            self.irrasjonal = S.One
            return
        if self.rasjonal == S.Zero or faktor.rasjonal == S.Zero: # x * 0 er 0
            self.rasjonal = S.Zero
            self.irrasjonal = S.One
            return

        resultat = Mul(self.rasjonal, self.irrasjonal, faktor.rasjonal, faktor.irrasjonal)
        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)
        
        
    def delt_med(self, dividend) -> None:
        """Deler self med dividend ved å gange self med 1/dividend

        Args:
            dividend (tall): Tallet som self deles med
        """
        if dividend.__class__ != tall: # Kan bare dele med et tall
            return
        a = dividend.kopi()
        a.resiprok()
        self.ganger(a)


    def fakultet(self) -> None:
        """Beregner fakultet av et naturlig tall: 1 * 2 * ... * self
        """        
        global varseltekst
        if self.irrasjonal != S.One or self.rasjonal < S.Zero or not self.rasjonal.is_integer:
            varseltekst = 'OBS! Kan bare beregne fakultet av naturlige tall.'
            return
        self.rasjonal = Integer(factorial(self.rasjonal))

  
    def fib(self) -> None:
        """Beregner fibonacci-tallet av et naturlig tall
        """        
        global varseltekst
        if self.irrasjonal != S.One or self.rasjonal < S.Zero or not self.rasjonal.is_integer:
            varseltekst = 'OBS! Kan bare beregne fibinacci-tallet av naturlige tall.'
            return
        self.rasjonal = fibonacci(self.rasjonal)

  
    def binom(self, k) -> None:
        """Beregner binominalkoeffisienten self over k

        Args:
            k (tall): nederste tall i binominalkoeffisienten n over k
        """        
        global varseltekst
        if k.__class__ != tall: # Kan bare regne med tall
            return
        if k.irrasjonal != S.One or not k.rasjonal.is_integer:
            varseltekst = 'OBS! Kan bare beregne binominalkoeffisienten over et heltall.'
            return
        if k.rasjonal < S.Zero:
            self.rasjonal = S.Zero
            self.irrasjonal = S.One
            return
        if self.irrasjonal == S.One:
            self.rasjonal = binomial(self.rasjonal, k.rasjonal)
        else:
            varseltekst = 'OBS! Kan bare beregne binominalkoeffisienten av rasjonale tall.'

  
    def opphøyd_i(self, potens) -> None:
        """Beregner self opphøyd i potens, dvs. opphøyer et tall i en potens

        Args:
            potens (tall): potensen som tallet opphøyes i
        """        
        global varseltekst
        if potens.__class__ != tall: # Kan bare regne med tall
            return
        if potens.rasjonal == S.Zero:
            self.rasjonal = S.One
            self.irrasjonal = S.One
            return
        if potens.rasjonal == S.One and potens.irrasjonal == S.One:
            return

        resultat = Pow(Mul(self.rasjonal, self.irrasjonal), Mul(potens.rasjonal, potens.irrasjonal))
        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)


    def lb(self) -> None:
        """Beregner logaritme til base 2
        """        
        global varseltekst
        if self.rasjonal < S.Zero:
            varseltekst = 'OBS! Kan bare beregne logaritmus av positive tall!'
            return
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne logaritmus av et tall!'
            return
        if self.rasjonal == S.Zero: # lb(0) = -∞
            self.rasjonal = S.NegativeInfinity
            self.irrasjonal = S.One
            return
        if self.rasjonal == S.Infinity: # lb(∞) = ∞
            return
        
        # Beregn logaritmen til base 2 ved hjelp av sympy
        resultat = log(Mul(self.rasjonal, self.irrasjonal), 2)
        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)
        

    def lg(self) -> None:
        """Beregner logaritme til base 10
        """        
        global varseltekst
        if self.rasjonal < S.Zero:
            varseltekst = 'OBS! Kan bare beregne logaritmus av positive tall!'
            return
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne logaritmus av et tall!'
            return
        if self.rasjonal == S.Zero: # lb(0) = -∞
            self.rasjonal = S.NegativeInfinity
            self.irrasjonal = S.One
            return
        if self.rasjonal == S.Infinity: # lb(∞) = ∞
            return
        
        # Beregn logaritmen til base 2 ved hjelp av sympy
        resultat = log(Mul(self.rasjonal, self.irrasjonal), 10)
        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)


    def ln(self) -> None:
        """Beregner logaritme til base e, den naturlige logaritme
        """        
        global varseltekst
        if self.rasjonal < S.Zero:
            varseltekst = 'OBS! Kan bare beregne logaritmus av positive tall!'
            return
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne logaritmus av et tall!'
            return
        if self.rasjonal == S.Zero: # lb(0) = -∞
            self.rasjonal = S.NegativeInfinity
            self.irrasjonal = S.One
            return
        if self.rasjonal == S.Infinity: # lb(∞) = ∞
            return
        
        # Beregn logaritmen til base 2 ved hjelp av sympy
        resultat = log(Mul(self.rasjonal, self.irrasjonal), 2)
        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)


    def sin(self) -> None:
        """Beregner sinus av en vinkel i grader eller radianer
        """        
        global varseltekst
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne sinus av et tall!'
            return
        if self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity:
            varseltekst = 'OBS! Kan ikke beregne sinus av uendelig!'
            self.rasjonal = S.NaN
            return
        
        # Beregn sinus ved hjelp av sympy
        if self.irrasjonal == S.Pi:
            resultat = sin(Mul(self.rasjonal, self.irrasjonal))
        else:
            resultat = sin(Mul(self.rasjonal, self.irrasjonal, S.Pi / 180))

        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)


    def cos(self) -> None:
        """Beregner kosinus av en vinkel i grader eller radianer
        """        
        global varseltekst
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne kosinus av et tall!'
            return
        if self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity:
            varseltekst = 'OBS! Kan ikke beregne kosinus av uendelig!'
            self.rasjonal = S.NaN
            return
        
        # Beregn kosinus ved hjelp av sympy
        if self.irrasjonal == S.Pi:
            resultat = cos(Mul(self.rasjonal, self.irrasjonal))
        else:
            resultat = cos(Mul(self.rasjonal, self.irrasjonal, S.Pi / 180))

        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)


    def tan(self) -> None:
        """Beregner tangens av en vinkel i grader eller radianer
        """        
        global varseltekst
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne tangens av et tall!'
            return
        if self.rasjonal == S.Infinity or self.rasjonal == S.NegativeInfinity:
            varseltekst = 'OBS! Kan ikke beregne tangens av uendelig!'
            self.rasjonal = S.NaN
            return
        
        # Beregn tangens ved hjelp av sympy
        if self.irrasjonal == S.Pi:
            if abs(self.rasjonal) % S.One == S.Half:
                varseltekst = 'OBS! Tangens av π/2 + nπ er ikke definert!'
                self.rasjonal = S.NaN
                self.irrasjonal = S.One
                return
            resultat = tan(Mul(self.rasjonal, self.irrasjonal))
        else:
            if abs(self.rasjonal) % Integer(180) == 90:
                varseltekst = 'OBS! Tangens av 90° + n*180° er ikke definert!'
                self.rasjonal = S.NaN
                self.irrasjonal = S.One
                return
            resultat = tan(Mul(self.rasjonal, self.irrasjonal, S.Pi / 180))

        self.rasjonal, self.irrasjonal = sjekk_resultat(resultat)


    def arcsin(self) -> None:
        """Beregner arcussinus av et tall i intervall [-1, 1]
        """        
        global varseltekst
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne arcsinus av et tall!'
            return
        if abs(Mul(self.rasjonal, self.irrasjonal)) > S.One:
            varseltekst = 'OBS! Tallet må være i intervall [-1, 1]!'
            self.rasjonal = S.NaN
            return
        
        # Beregn arcsinus ved hjelp av sympy
        resultat = asin(Mul(self.rasjonal, self.irrasjonal))

        self.rasjonal, self.irrasjonal = sjekk_resultat(Mul(resultat, 180 / S.Pi))


    def arccos(self) -> None:
        """Beregner arcuscosinus av et tall i intervall [-1, 1]
        """        
        global varseltekst
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne arccosinus av et tall!'
            return
        if abs(Mul(self.rasjonal, self.irrasjonal)) > S.One:
            varseltekst = 'OBS! Tallet må være i intervall [-1, 1]!'
            self.rasjonal = S.NaN
            return
        
        # Beregn arccosinus ved hjelp av sympy
        resultat = acos(Mul(self.rasjonal, self.irrasjonal))

        self.rasjonal, self.irrasjonal = sjekk_resultat(Mul(resultat, 180 / S.Pi))


    def arctan(self) -> None:
        """Beregner arctangens av et tall
        """        
        global varseltekst
        if self.rasjonal == S.NaN:
            varseltekst = 'OBS! Kan bare beregne tangens av et tall!'
            return
        
        # Beregn arctangens ved hjelp av sympy
        resultat = atan(Mul(self.rasjonal, self.irrasjonal))

        self.rasjonal, self.irrasjonal = sjekk_resultat(Mul(resultat, 180 / S.Pi))


    def utvid(self) -> None:
        """Utvider et uttrykk ved å bruke sympys expand(func=True)
        """
        if self.rasjonal != S.NaN and self.irrasjonal != S.One:
            uttrykk = expand(Mul(self.rasjonal, self.irrasjonal), func=True)
            self.rasjonal, self.irrasjonal = sjekk_resultat(uttrykk)


    def evaluer(self) -> None:
        """Beregner et irrasjonalt tall til et rasjonalt tall ved å bruke sympy
        """
        if self.rasjonal != S.NaN and self.irrasjonal != S.One:
            tall = Mul(self.rasjonal, self.irrasjonal)
            self.rasjonal = Rational(tall.evalf(n=MAKS_SIFRE, chop=Pow(10, -MAKS_SIFRE)))
            self.irrasjonal = S.One
            self.feil = Integer(log(abs(Add(tall, -self.rasjonal)), 10).evalf()) # Feil i beregningen
            if self.feil < -100:
                self.feil = S.Zero


def beregne(tallene: list[tall], operasjon: str) -> None:
    """
    Anvender operasjonen på tallene i stabelen. Resultatet legges på stabelen

    Args:
        tallene (list[tall]): stabelen som inneholder noen tall
        operasjon (str): operasjonen som skal anvendes
    """
    global varseltekst
    if len(tallene) == 0:
        return
    match operasjon:
        case 'slett':
            tallene.pop()
        case 'bytt':
            if len(tallene) > 1:
                a = tallene[-2]
                tallene[-2] = tallene[-1]
                tallene[-1] = a
        case 'kopi':
            tallene.append(tallene[-1].kopi())
        case '+':
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].pluss(a)
        case '-':
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].minus(a)
            elif len(tallene) == 1:
                tallene[-1].rasjonal *= S.NegativeOne # Bytter fortegn
        case '--': # bytt fortegn
            tallene[-1].rasjonal *= S.NegativeOne # Bytter fortegn
        case '*':
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].ganger(a)
        case '/':
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].delt_med(a)
        case '%': # Det er prosent (ikke modulo-operasjonen)
            if len(tallene) > 1:
                a = tallene.pop()
                a.delt_med(tall('100'))
                tallene[-1].ganger(a)
        case '!': # fakultet
            tallene[-1].fakultet()
        case '()': # binominalkoeffisient
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].binom(a)
        case 'binom': # binominalkoeffisient
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].binom(a)
        case '**': # opphøye et tall
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].opphøyd_i(a)
        case 'v': # kvadratrot
            tallene[-1].opphøyd_i(tall('1/2'))
        case 'rot': # a-te rot
            if len(tallene) > 1:
                a = tallene.pop()
                a.resiprok()
                tallene[-1].opphøyd_i(a)
        case 'lg': # logaritme til base 10
            tallene[-1].lg()
        case 'ln': # naturlig logaritme (base e)
            tallene[-1].ln()
        case 'lb': # logaritme til base 2
            tallene[-1].lb()
        case 'sin':
            tallene[-1].sin()
        case 'cos':
            tallene[-1].cos()
        case 'tan':
            tallene[-1].tan()
        case 'arcsin':
            tallene[-1].arcsin()
        case 'arccos':
            tallene[-1].arccos()
        case 'arctan':
            tallene[-1].arctan()
        case 'fib': # bergne n-te fibinacci-tall
            tallene[-1].fib()
        case '//': # Det er 1/x, ikke Pythons //
            tallene[-1].resiprok()
        case 'resiprok':
            tallene[-1].resiprok()
        case 'mod':
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].mod(a)
        case 'rest':
            if len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].mod(a)
        case 'hel':
            tallene[-1].gjør_hel()
        case 'grad': # omgjør fra radianer til grader
            a = tall('180')
            a.delt_med(tall(PI))
            tallene[-1].ganger(a)
        case 'rad': # omgjør fra grader til radianer
            a = tall(PI)
            a.delt_med(tall('180'))
            tallene[-1].ganger(a)
        case 'utvid': # utvider et uttrykk
            if tallene[-1].irrasjonal != S.One:
                tallene[-1].utvid()
        case 'eval': # evaluerer et uttrykk
            if tallene[-1].irrasjonal != S.One:
                tallene[-1].evaluer()
        case _:
            varseltekst = 'OBS! Ukjent operasjon'


def notasjon(tallet: Rational, feil: Integer) -> str:
    """Konverterer en brøk til vitenskapelig notasjon

    Args:
        tallet (Rational): tallet som skal skrives i vitenskapelig notasjon
        feil (Integer): feil i beregningen
    Returns:
        str: tallet i vitenskapelig notasjon slik at den kan brukes i LaTeX
    """
    if not isinstance(tallet, Rational):
        return ''
    fortegn: str = '-' if tallet < S.Zero else ''
    ellipsis: str = ''
    teller: Integer = abs(numer(tallet))
    nevner: Integer = denom(tallet)
    scientific_str: str = ''

    if teller < nevner:
        # Beregn hvor mye mindre enn 1 tallet er, dvs. hvor mange nuller bak kommaet det har
        rest: Integer = teller
        nuller: Integer = S.Zero

        # Tell hvor ofte vi må gange rest med ti til det blir større enn nevner
        while rest < nevner:
            rest *= Integer(10)
            nuller += S.One
        
        # Gjør en lang divisjon bak nullene
        huskelapp: dict = {}  # rest -> indeks i desimaldelen
        desimaldel: list[str] = []
        indeks: int = 0
        # Beregn sifrene til vi finner perioden eller når ønsket presisjon
        while (rest not in huskelapp) and (len(desimaldel) < MAKS_SIFRE):
            siffer = rest // nevner
            desimaldel.append(str(siffer))
            huskelapp[rest] = indeks # Husk at vi fant rest ved indeks
            rest %= nevner # Beregne ny rest
            rest *= Integer(10) 
            indeks += S.One
        
        if len(desimaldel) >= MAKS_SIFRE: # Vi har nådd ønsket presisjon, men ikke funnet en periode
            førperiode = ''.join(desimaldel)
            periode = ''
            feil = max(-MAKS_SIFRE, feil)
            ellipsis = '\\ldots'
        else: # Vi har funnet perioden fordi resten står på huskelappen
            periode_start_indeks = huskelapp[rest]
            førperiode = ''.join(desimaldel[:periode_start_indeks])
            periode = ''.join(desimaldel[periode_start_indeks:])
                
        # Sjekk om der er 0- eller 9-periode
        if periode == '0': # 0-periode betyr at tallet er nøyaktig
            periode = ''
        if periode == '9': # 0,9999... = 1
            # Vi må øke siste sifre før perioden med 1. Derfor må vi omgjøre det til int
            sifrene = str(int('1' + førperiode) + 1) # For å ungå problemer med ledende nuller, tilføyer vi '1' på starten
            if sifrene[0] == '1':
                førperiode = sifrene[1:] # ta bort '1' som vi har tilføyd
            else: # Det var en i mente slik at '1' blir '2'
                førperiode = '1' + sifrene[1:]
                nuller -= S.One # alt flyttes en plass til venstre
            periode = ''

        if nuller == 1: # Det er bare en null foran kommaet
            if periode:
                scientific_str = '0,' + førperiode + '\\overline{' + periode + '}'
            elif førperiode:
                scientific_str = '0,' + førperiode + ellipsis
            else:
                scientific_str = '0'
        elif nuller < MAKS_SIFRE//Integer(2): # Hvis det er ikke altfor mange nuller, kan tallet framstilles som desimaltall
            if periode:
                if førperiode:
                    scientific_str = '0,' + ('0' * int(nuller - 1)) + førperiode + '\\overline{' + periode + '}'
                    # OBS: Det er (nuller - 1) fordi en null står foran kommaet
                elif periode[-1] == '0': # Flytt nuller fra slutten av perioden til begynnelsen
                    uten_nuller = periode.rstrip('0') # Det kan være flere 0er
                    n_nuller = len(periode) - len(uten_nuller) # Beregn hvor mange nuller det er
                    scientific_str = '0,' + ('0' * int(nuller - n_nuller - S.One)) + '\\overline{' + ('0' * n_nuller) + uten_nuller + '}'
                else: # Det er ikke nuller som må flyttes
                    scientific_str = '0,' + ('0' * int(nuller - 1)) + '\\overline{' + periode + '}'
            elif førperiode:
                scientific_str = '0,' + ('0' * int(nuller - 1)) + førperiode + ellipsis
            else:
                scientific_str = '0'
        else: # Tallet må framstilles i vitenskapelig notasjon
            if feil != S.Zero:
                feil -= nuller
            if førperiode: # Det er sifre foran perioden
                første_siffer = førperiode[0]
                førperiode = førperiode[1:]
                if periode:
                    scientific_str = første_siffer + ',' + førperiode + '\\overline{' + periode + '}'
                elif førperiode:
                    scientific_str = første_siffer + ',' + førperiode + ellipsis + ' \\cdot 10^{-' + str(nuller) + '}'
                elif første_siffer == '1':
                    scientific_str = '10^{-' + str(nuller) + '}'
                else:
                    scientific_str = første_siffer + ' \\cdot 10^{-' + str(nuller) + '}'
            elif periode: # Det er ikke sifre foran perioden
                første_siffer = periode[0]
                if len(periode) > 1:
                    periode = periode[1:] + første_siffer
                scientific_str = første_siffer + ',\\overline{' + periode + '} \\cdot 10^{-' + str(nuller) + '}'
    else: # teller >= nevner, dvs, tallet >= 1
        # Beregn først tallet foran kommaet
        heldel = Integer(abs(tallet))
        helstr = str(heldel)
        eksponent = len(helstr) - 1  # eksponenten er antall sifre i heldelen minus 1 (minus 1 fordi det skal være ett siffer foran kommaet)
        # Beregn desimaldelen ved lang divisjon
        rest = teller % nevner * Integer(10)
        huskelapp: dict = {}  # rest -> indeks i desimaldelen
        desimaldel: list[str] = []
        indeks: int = 0
        # Beregn sifrene til vi finner perioden eller når ønsket presisjon
        while (rest not in huskelapp) and (len(desimaldel) < MAKS_SIFRE):
            siffer = rest // nevner
            desimaldel.append(str(siffer))
            huskelapp[rest] = indeks # Husk at vi fant rest ved indeks
            rest %= nevner # Beregne ny rest
            rest *= Integer(10)
            indeks += S.One
            
        if len(desimaldel) >= MAKS_SIFRE: # Vi har nådd ønsket presisjon, men ikke funnet en periode
            førperiode = ''.join(desimaldel)
            periode = ''
            feil = max(-MAKS_SIFRE, feil)
            førperiode = førperiode + '\\ldots'
        else: # Vi har funnet perioden fordi resten står på huskelappen
            periode_start_indeks = huskelapp[rest]
            førperiode = ''.join(desimaldel[:periode_start_indeks])
            periode = ''.join(desimaldel[periode_start_indeks:])
                    
        # Sjekk om der er 0- eller 9-periode
        if periode == '0': # 0-periode betyr at tallet er nøyaktig
            periode = ''
        if periode == '9': # 0,9999... = 1
            # Vi må øke siste sifre før perioden med 1. Derfor må vi omgjøre det til int
            sifrene = str(int('1' + førperiode) + 1) # For å ungå problemer med ledende nuller, tilføyer vi '1' på starten
            if sifrene[0] == '1':
                førperiode = sifrene[1:] # ta bort '1' som vi har tilføyd
            else: # Det var en i mente slik at '1' blir '2'
                førperiode = '1' + sifrene[1:]
                eksponent -= S.One # alt flyttes en plass til venstre
            periode = ''

        if eksponent <= S.One: # Hvis heldelen har bare en eller to sifre, kan de bli stående
            første_siffer = helstr
            brøkdel = førperiode + '\\overline{' + periode + '}' if periode else førperiode
            if brøkdel:
                scientific_str = første_siffer + ',' + brøkdel
            else:
                scientific_str = første_siffer
        elif eksponent < MAKS_SIFRE//Integer(4): # Hvis tallet ikke er for stor, kan den framstilles som desimaltall
            første_siffer = '{:,}'.format(int(heldel)).replace(',', '~')
            brøkdel = førperiode + '\\overline{' + periode + '}' if periode else førperiode
            eksponent = S.Zero
            if brøkdel:
                scientific_str: str = første_siffer + ',' + brøkdel
            else:
                scientific_str: str = første_siffer
        else: # ellers bruker vi vitenskapelig notasjon med ett siffer foran kommaet
            if feil != S.Zero:
                feil += eksponent
            første_siffer = helstr[0]
            førperiode = helstr[1:] + førperiode
            if periode:
                førperiode = førperiode.replace(periode, 'x').rstrip('x')
                while førperiode and førperiode[-1] == periode[-1]:
                    periode = førperiode[-1] + periode[:-1]
                    førperiode = førperiode[:-1]
                brøkdel = førperiode + '\\overline{' + periode + '}'
            else:
                brøkdel = førperiode
            if brøkdel:
                scientific_str = første_siffer + ',' + brøkdel + ' \\cdot 10^{' + str(eksponent) + '}'
            elif første_siffer == '1':
                scientific_str = '10^{' + str(eksponent) + '}'
            else:
                scientific_str = første_siffer + ' \\cdot 10^{' + str(eksponent) + '}'

    if feil != S.Zero:
        return '\\approx ' + fortegn + scientific_str + '\\qquad \\pm10^{' + str(feil) + '}'
    return '= ' + fortegn + scientific_str 


def heltall_til_latex(verdi: Integer, feil: Integer) -> str:
    """Konverterer et heltall til en LaTeX-streng.
    Args:
        verdi (Integer): Heltallet som skal konverteres til LaTeX-format.
        feil (Integer): Feil i beregningen
    Returns:
        str: LaTeX-strengen som representerer tallet
    """
    likhetstegn: str = '= '
    ellipsis: str = ''
    feilstr: str = '' if feil == S.Zero else '\\qquad \\pm10^{' + str(feil) + '}'

    eksponent = num_digits(verdi) - 1

    sifre = digits(verdi) # Tallets sifre som liste
    første_siffer = str(-sifre[1]) if sifre[0] < 0 else str(sifre[1])  # Første siffer i tallet

    # Resten av sifrene i tallet, maks MAKS_SIFRE sifre etter første siffer
    if len(sifre) > MAKS_SIFRE + 2: # Det er flere sifre enn det vi kan vise
        resten_av_tall = ''.join(str(sifr) for sifr in sifre[2:MAKS_SIFRE + 2])
        ikke_bare_nuller = next((i for i, x in enumerate(sifre[MAKS_SIFRE + 2:]) if x != 0), -1)
        if ikke_bare_nuller != -1: # Det finnes sifre som ikke vises. Derfor er resultatet ikke eksakt.
            likhetstegn = '\\approx '
            ellipsis = '\\ldots'
            feileksponent = max(eksponent - MAKS_SIFRE - ikke_bare_nuller - 1, feil)
            feilstr = '\\qquad \\pm10^{' + str(feileksponent) + '}'
    else:
        resten_av_tall = ''.join(str(sifr) for sifr in sifre[2:]) 

    uten_nuller = resten_av_tall.rstrip('0') # Fjerner nuller på slutten av resten av tallet

    komma = ',' if len(uten_nuller) > 0 else '' # Hvis det er ingen sifre etter kommaet, fjerner vi kommaet

    return likhetstegn + første_siffer + komma + uten_nuller + ellipsis + '\\cdot 10^{' + str(eksponent) + '}' + feilstr


def til_latex(tallet: tall) -> str:
    """Konverterer et tall til en LaTeX-streng.
    Args:
        tallet (tall): Tallet som skal konverteres til LaTeX-format

    Returns:
        str: LaTeX-strengen som representerer tallet
    """
    if tallet.rasjonal == S.Infinity:
        return '= \\infty'
    if tallet.rasjonal == S.NegativeInfinity:
        return '= -\\infty'
    if tallet.rasjonal == S.NaN:
        return '\\bot'

    verdi = tallet.rasjonal if isinstance(tallet.rasjonal, Rational) else None

    limit = Pow(10, MAKS_SIFRE)
    if verdi and tallet.irrasjonal == S.One: # Rasjonale tall kan vises som som brøk, blandet tall eller heltall
        if verdi.is_integer: # Det er et heltall
            if abs(verdi) < limit:
                # Formater tallet slik at de er lettere å lese, dvs. 1000000 blir 1 000 000
                if tallet.feil != 0:
                    return '\\approx {:,}'.format(int(verdi)).replace(',', '~') + '\\qquad \\pm10^{' + str(tallet.feil) + '}'
                else:
                    return '= {:,}'.format(int(verdi)).replace(',', '~')
            # Et heltall med flere enn MAKS_SIFRE sifre vises i vitenskapelig notasjon
            return heltall_til_latex(Integer(verdi), tallet.feil)
        if verdi.is_rational: # Det er en brøk
            teller: Integer = numer(verdi) # Henter teller fra brøken
            nevner: Integer = denom(verdi) # Henter nevner fra brøken
            if verdi > 1 or verdi < -1:
                if max(teller, nevner) < limit: # Det er et blandet tall som kan vises på vanlig måte og vitenskapelig notasjon
                    hel = Integer(verdi)
                    # Formater hel, teller og nevner slik at de er lettere å lese, dvs. 1000000 blir 1 000 000
                    nevnerstr = '{:,}'.format(int(nevner)).replace(',', '~')
                    tellerstr = '{:,}'.format(int(abs(Add(teller, -Mul(hel, nevner))))).replace(',', '~')
                    helstr = '= {:,}'.format(int(hel)).replace(',', '~')
                    return helstr + ' \\frac{' + tellerstr + '}{' + nevnerstr + '} ' + notasjon(verdi, tallet.feil)
                return notasjon(verdi, tallet.feil) 
            teller = abs(teller)
            if max(teller, nevner) < limit: # Det er en brøk som kan vises på vanlig måte
                # Formater teller og nevner slik at de er lettere å lese, dvs. 1000000 blir 1 000 000
                tellerstr = '{:,}'.format(int(teller)).replace(',', '~')
                nevnerstr = '{:,}'.format(int(nevner)).replace(',', '~')
                if verdi < 0: # Hvis tallet er negativt, må vi vise det med minus foran
                    return '= -\\frac{' + tellerstr + '}{' + nevnerstr + '}' + notasjon(verdi, tallet.feil)
                else:
                    return '= \\frac{' + tellerstr + '}{' + nevnerstr + '}' + notasjon(verdi, tallet.feil)
            # Hvis teller eller nevner er for store, brukes bare vitenskapelig notasjon
            return notasjon(verdi, tallet.feil)
        return '\\bot' # Det er ikke et tall som kan formateres til LaTeX-formatet

    # Tallet er irrasjonalt
    return '= ' + latex(Mul(verdi, tallet.irrasjonal), decimal_separator='comma', max=MAKS_SIFRE)


def skriv_resultat(tallet: tall) -> PIL.Image:
    """Lager et bilde som inneholder resultatet i LaTeX-formatet og viser det i vinduet.

    Args:
        tallet (tall): Tallet som er resultatet som skal vises

    Returns:
        ImageFile: Bildet som skal vises
    """
    # Lage en matplotlib-figur
    fig = Figure(figsize=(12, 0.8), dpi=100)
    ax = fig.add_subplot(111)
    ax.axis('off')
        
    # Lage en LaTeX-representasjon av det matematiske uttrykket
    uttrykk = '$' + til_latex(tallet) + '$'

       
    # Tilføy uttrykket som tekst til figuren
    ax.text(0.5, 0.5, uttrykk, fontsize=16, ha='center', va='center')

    # Trykke figuren på lerret
    lerret = FigureCanvasAgg(fig)
    buffer = io.BytesIO()
    lerret.print_png(buffer)
    buffer.seek(0)
    bildet = PIL.Image.open(buffer)
    return bildet

        
def stakke(stabel: list[tall]) -> None:
    """Oppdaterer visningene av tallene i stabelen

    Args:
        stabel (list[tall]): Liste av tallene som skal vises
    """    
    # Oppdatere skjerm-stabelen med PIL-bildene
    for i in range(8):
        if i < len(stabel):
            # I stabel står nederste tall bakest, men i lbl_stabel har nederste linje indeks 7
            img = skriv_resultat(stabel[- (i + 1)])
            lbl_stabel[7 - i].image = toga.Image(img)
        else:
            lbl_stabel[7 - i].image = None
    resultat = stabel[-1].kopi()
    resultat.evaluer()
    img = skriv_resultat(resultat)
    lbl_resultat.image = toga.Image(img)
    if resultat.feil == S.Zero:
        lbl_stabel[7].image = toga.Image(img)


class SymKalkulator(toga.App):
    """Det er klassen som lager app'en 'SymKalkulator'.

    Args:
        toga (App): Toga er GUI'en som brukes
    """    
    def startup(self) -> None:
        """Startup er en funksjon som hver Toga-app må ha. Det er funksjonen som kjøres når SymKalkulator() kalles.
        """        
        global varseltekst
        stabel: list[tall] = []

        async def inntastet(self) -> None:
            """Det er en funksjon som kjører asynkron i bakgrunnen.
               Den sjekker bestandig om noe har blitt tastet inn i feltet nederst i app'en.
               Når noe har blitt tastet inn, analyserer den hva som har blitt tastet inn og avgjør den hva som skal gjøres med det.
            """            
            global varseltekst
            innput: str = self.value # self.value er det som har blitt tastet inn som string
            self.value = '' # Etter verdien har blitt lagret i variabelen 'innput' tilbakestilles feltet og venter på nye ting som kan tastes inn
            if len(innput) == 0: # Hvis det var ikke noe som har blitt tastet inn, skjer ingenting. App'en venter til noe blir tastet inn
                return
            if innput[-1].isdigit(): # Ett tall har blitt tastet inn
                nytall = tall(innput)
                stabel.append(nytall) # Tallet legges på stabelen
            elif len(innput) == 1: # Ikke ett tall, men bare ett tegn har blitt tastet inn
                if innput == 'e': # Det kan være Eulers tall e
                    nytall = tall('e')
                    stabel.append(nytall) # Tallet legges på stabelen
                else: # Eller det er en operasjon fra listen ['+', '-', '*', '/', '%', '!', 'v']
                    operasjon = innput
                    beregne(stabel, operasjon) # Anvender operasjonen på tallene som ligger på stabelen
            elif innput in {'-e', 'pi', '-pi', 'tau', '-tau', 'fi', '-fi', 'phi', '-phi', 'oo', '-oo'}: # Det har blitt tastet inn flere enn ett tegn
                nytall = tall(innput)
                stabel.append(nytall)
            elif (innput[-1] in ['+', '-', '*', '/', '%', '!', 'v']) and innput[-2].isdigit():
                # Det er ett tall med en operasjon som består av ett tegn
                operasjon = innput[-1]
                nytall = tall(innput[:-1])
                stabel.append(nytall)
                beregne(stabel, operasjon)
            elif len(innput) > 2 and innput[-2:] == '**' and innput[-3].isdigit():
                # Det er ett tall med operasjonen ** (opphøyd i)
                operasjon = '**'
                nytall = tall(innput[:-2])
                stabel.append(nytall)
                beregne(stabel, operasjon)
            else:
                operasjon = innput
                beregne(stabel, operasjon)
            stakke(stabel)
            if varseltekst != '': # Det har oppstått en feil under beregningen, f.eks. at operasjonen er ukjent
                varseldialog = toga.ErrorDialog('Varsel', varseltekst)
                await SymKalkulator.dialog(self, varseldialog)
                varseltekst = ''


        async def kopier(self) -> None:
            """Det er en asynkron funksjon som kalles når kopier-funskjonen brukes fra Commands-menyen eller ved å trykke Ctrl+C
            """            
            if len(stabel) > 0: # Den kopier siste tall fra stabelen hvis det er et tall på stabelen
                tallet = stabel[-1].kopi()
                pyperclip.copy(str(Mul(tallet.rasjonal, tallet.irrasjonal).evalf(n=MAKS_SIFRE)))

        # Det følgende er fortsatt del av startup-funksjonen. Her defineres hvordan vinduet ser ut og fungerer.
        innhold = toga.Box(style=Pack(direction=COLUMN, padding=10)) # Innholdet i hovedvinduet
        self.main_window = toga.MainWindow(title=self.formal_name) # Hovedvinduet må hete 'main_window'
        self.main_window.content = innhold
        # Jeg bruker norske betegnelser i menyen
        self.commands[toga.Command.VISIT_HOMEPAGE].text = 'Besøk hjemmesida'
        self.commands[toga.Command.ABOUT].text = 'Om SymKalkulator'
        self.commands[toga.Command.EXIT].text = 'Avslutt'
        hjelp_cmd = toga.Command(
            self.hjelp,
            text='Hjelp',
            tooltip='Vis funksjonene som kan brukes',
            shortcut=toga.Key.MOD_1 + 'h',
            group=toga.Group.HELP
        )
        copy_cmd = toga.Command(
            kopier,
            text='Kopier',
            tooltip='Kopier resultat til clipboard',
            shortcut=toga.Key.MOD_1 + 'c'
        )
        self.commands.add(hjelp_cmd)
        self.commands.add(copy_cmd)
        
        # Lage stabel-boxen
        stabel_box = toga.Box(style=Pack(direction=COLUMN, height=665, padding=2))
        global lbl_stabel
        lbl_stabel = [toga.ImageView(None, style=Pack(height=80, padding=2)) for _ in range(8)]
        for i in range(7):
            stabel_box.add(lbl_stabel[i])
        stabel_box.add(toga.Divider())
        stabel_box.add(lbl_stabel[7])

        # Lage innput-boxen
        input_box = toga.Box(style=Pack(direction=COLUMN, alignment=CENTER, padding=2))
        self.innput = toga.TextInput(style=Pack(text_align=CENTER, width=1200, height=50, padding=2, font_size=20), on_confirm=inntastet)
        input_box.add(self.innput)

        # Lage resultat-boxen
        resultat_box = toga.Box(style=Pack(direction=COLUMN, alignment=CENTER, padding=2))
        global lbl_resultat
        lbl_resultat = toga.ImageView(None, style=Pack(height=80, padding=2))
        resultat_box.add(lbl_resultat)

        # Tilføy alle deler til hovedboxen
        innhold.add(stabel_box)
        innhold.add(toga.Divider())
        innhold.add(toga.Divider())
        innhold.add(input_box)
        innhold.add(toga.Divider())
        innhold.add(resultat_box)
        innhold.add(toga.Divider())
        innhold.add(toga.Divider())
        self.main_window.show()

    async def hjelp(self, widget):
        hjelpedialog = toga.StackTraceDialog('Hjelpevindu', hjelpetittel, hjelpetekst)
        await self.main_window.dialog(hjelpedialog)


def main():
    return SymKalkulator()
