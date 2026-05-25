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
import sympy as sp
from sympy.ntheory.digits import digits
import sympy.physics.units as u
import sympy.physics.units.util as uu
from sympy.parsing.latex import parse_latex

import re # Regular Expressions for å analysere det som blir tastet inn
import pyperclip # for å kopiere verdier til utklippstavla

MAKS_SIFRE = sp.Integer(32) # Maksimalt antall sifre som kan vises i resultatfeltet

# Stringsymboler
UENDELIG = u'\u221E'
PI = u'\u03C0'
TAU = u'\u03C4'
FI = u'\u03C6'

# Sympy-atomer
ETTALLET = sp.S.One
NULLTALLET = sp.S.Zero
IKKETALL = sp.S.NaN

# Lage hjelpetekst

kommandoer =   {'rydd':     'Rydd stabelen',
                'slett':    'Slett y fra stabelen',
                'bytt':     'Bytt x og y',
                'kopi':     'Kopier y',
                '+':        'Adder x + y',
                '-':        'Subtraher x - y',
                '--':       'Bytt fortegn av y',
                '*':        'Multipliser x·y',
                '/':        'Divider x/y',
                '%':        'Beregn y prosent av x',
                '!':        'Beregn fakultet av y',
                '()':       'Beregn binominalkoeffisient',
                '**':       'Beregn potens xʸ',
                'v':        'Beregn kvadratrot av y',
                'rot':      'Beregn y-te rot av x',
                'lg':       'Beregn briggsk logaritme av y',
                'ln':       'Beregn naturlig logaritme av y',
                'lb':       'Beregn binær logaritme y',
                'sin':      'Beregn sinus av y',
                'cos':      'Beregn cosinus av y',
                'tan':      'Beregn tangens av y',
                'arcsin':   'Beregn arcsinus av y',
                'arccos':   'Beregn arccosinus av y',
                'arctan':   'Beregn arctangens av y',
                'fib':      'Beregn det y-te fibonacci-tallet',
                '//':       'Beregn 1/y',
                'resiprok': 'Beregn 1/y',
                'mod':      'Beregn rest av divisjonen x/y',
                'rest':     'Beregn rest av divisjonen x/y',
                'sum':      'Beregn summen av alle tallene i stabelen',
                'Ø':        'Beregn gjennomsnitt av alle tallene i stabelen',
                'prod':     'Beregn produktet av alle tallene i stabelen',
                'hel':      'Omgjør y til heltall',
                '<enhet>':  'Omgjør y til enheten <enhet>',
                'SI':       'Omgjør y til SI-enheter',
                'utvid':    'Utvid et uttrykk',
                'eval':     'Evaluer et uttrykk'} # Det er mange kommandoer som kan brukes i kalkulatoren, og de har alle en beskrivelse som vises i hjelpeteksten

spesielle_tall = {'e': 'Eulers tall e',
                  'π': 'sirkeltallet π',
                  'pi': 'sirkeltallet π',
                  'tau': 'Det dobbelte sirkeltallet ' + TAU,
                  'fi': 'Det irrasjonale tallet ' + FI,
                  'phi': 'Det irrasjonale tallet ' + FI,
                  'oo': 'uendelig ' + UENDELIG,
                  'G': 'gravitasjonskonstanten G',
                  'c': 'lysets hastighet c',
                  'hbar': 'Plancks redusete konstant ħ',
                  'ħ': 'Plancks redusete konstant ħ',
                  'k': 'Coulombs konstant k' + u'\u2091',
                  'R': 'gasskonstanten R'} # Det er spesielle tall som kan tastes inn uten å være en del av en operasjon, f.eks. 'e' eller 'pi'

hjelpetittel = 'Hvis x og y er nederste tallene, for <funksjon>       tast <innput>\n\n'
hjelpetekst = ''
for kommando in kommandoer:
    hjelpetekst += f'{kommandoer[kommando]:<40} \u2192 {kommando:^6}\n'
for tall in spesielle_tall:
    hjelpetekst += f'{spesielle_tall[tall]:<40} \u2192 {tall:^6}\n'

varseltekst = '' # Global variable for varsler


def split_måltall_enhet(enhet: sp.Mul) -> tuple[sp.Rational, sp.Expr, sp.Mul]:
    """Denne funksjonen tar en enhet og deler den opp i måltall og enhet ved å bruke sympy. Den forenkler også enheten ved å bruke sympy.
    
    Args:
        enhet (sp.Mul): Enheten som skal deles opp, f.eks. 2*u.meter/u.second
    
    Returns:
        sp.Rational:   Måltallet i enheten, f.eks. 2 i 2*u.meter/u.second
        sp.Expr:       Den irrasjonale delen av måltallet, f.eks. 1 i 2*u.meter/u.second
        sp.Mul:        Enheten uten måltallet, f.eks. u.meter/u.second i 2*u.meter/u.second
    """
    uttrykk = uu.quantity_simplify(enhet)
    rasjonale_faktorer = []
    irrasjonale_faktorer = []
    enhet_faktorer = []

    for faktor in sp.Mul.make_args(uttrykk):
        if faktor.has(u.Quantity):
            enhet_faktorer.append(faktor)
        elif isinstance(faktor, sp.Rational):
            rasjonale_faktorer.append(faktor)
        elif isinstance(faktor, sp.Float):
            rasjonale_faktorer.append(sp.Rational(str(faktor)))
        else:
            irrasjonale_faktorer.append(faktor)

    rasjonal = sp.Mul(*rasjonale_faktorer) if rasjonale_faktorer else ETTALLET
    irrasjonal = sp.Mul(*irrasjonale_faktorer) if irrasjonale_faktorer else ETTALLET
    enhet = sp.Mul(*enhet_faktorer) if enhet_faktorer else ETTALLET

    return rasjonal, irrasjonal, enhet


def samme_dimensjon(enhet1: sp.Mul, enhet2: sp.Mul) -> bool:
    """Denne funksjonen sammenligner to enheter og sjekker om de har samme dimensjon.

    Args:
        enhet1 (sp.Mul): Den første enheten som skal sammenlignes.
        enhet2 (sp.Mul): Den andre enheten som skal sammenlignes.

    Returns:
        bool: True hvis enhetene er like, False ellers.
    """
    _, _, forenklet_enhet1 = split_måltall_enhet(enhet1)
    _, _, forenklet_enhet2 = split_måltall_enhet(enhet2)
    enhet2 = u.convert_to(forenklet_enhet2, forenklet_enhet1)
    _, _, forenklet_enhet2 = split_måltall_enhet(enhet2)
    return forenklet_enhet1 == forenklet_enhet2


def parse_enhet(tekst: str) -> sp.Mul:
    """Denne funksjonen tar en string og prøver å parse den som en sympy-enhet. Hvis det ikke er en gyldig enhet, returneres ETTALLET.

    Args:
        tekst (str): Teksten som skal parses som en enhet, f.eks. 'm/s' for meter per sekund

    Returns:
        sp.Mul: Enheten som teksten representerer, f.eks. u.meter/u.second for 'm/s'
    """
    global varseltekst

    prefixes = {'T': 1e12, 'G': 1e9, 'M': 1e6, 'k': 1e3, 'h': 1e2, 'd': 1e-1, 'c': 1e-2, 'm': 1e-3, 'µ': 1e-6, 'n': 1e-9, 'p': 1e-12}
    enheter = {'1': ETTALLET, 'min': u.minute, "'": u.minute, '"': u.second, 'a': u.year, 'd': u.day, 'grad': u.degree, 'in': u.inch, 'Å': u.angstrom} # Det er noen enheter som kan tastes inn uten å være en del av en operasjon, f.eks. 'min' for minutter eller 'in' for tommer
    
    if tekst == '': # Det er ingen enhet
        varseltekst = 'OBS! Ugyldig innput. Det er ingen enhet.'
        return ETTALLET
    if tekst in enheter: # Sjekker om teksten er en gyldig enhet
        return enheter[tekst]
    if tekst in u.__all__[62:358]: # Sjekker om teksten er en gyldig sympy-enhet
        return eval('u.' + tekst)
    tekst = tekst.replace(' ', '') # Fjerner mellomrom fra teksten
    tekst = tekst.replace('**', '^') # erstatter ** med ^ for å gjøre det lettere å parse potens av enheter
    if tekst[0] == '(' and tekst[-1] == ')': # Sjekker om teksten er et uttrykk i parentes, f.eks. '(meter/second)**2'
        return parse_enhet(tekst[1:-1])
    if '/' in tekst: # Sjekker om teksten er en brøk av enheter, f.eks. 'meter/second'
        teller, nevner = tekst.split('/')
        return parse_enhet(teller) / parse_enhet(nevner)
    if '*' in tekst: # Sjekker om teksten er et produkt av enheter, f.eks. 'newton*meter'
        faktorer = tekst.split('*')
        resultat = ETTALLET
        for faktor in faktorer:
            resultat *= parse_enhet(faktor)
        return resultat
    if '^' in tekst: # Sjekker om teksten er en potens av enheter, f.eks. 'meter**2'
        base, eksponent = tekst.split('^')
        return parse_enhet(base)**sp.Rational(eksponent)
    if tekst[0] in prefixes: # Sjekker om teksten er en gyldig SI-prefiks, f.eks. 'k' for kilo
        prefix = tekst[0]
        enhet = tekst[1:]
        return prefixes[prefix] * parse_enhet(enhet)
    
    varseltekst = 'OBS! Ugyldig innput. Enheten er ukjent.'
    return ETTALLET
    

def latex_enhet(enhet: sp.Mul) -> str:
    """Denne funksjonen tar en sympy-enhet og returnerer en LaTeX-representasjon av den.

    Args:
        enhet (sp.Mul): Enheten som skal konverteres til LaTeX, f.eks. u.meter/u.second

    Returns:
        str: LaTeX-representasjonen av enheten, f.eks. '\\frac{m}{s}' for u.meter/u.second
    """
    abbrevs = {u.hour: '\\text{h}', u.minute: '\\text{min}', u.year: '\\text{a}', u.day: '\\text{d}', 
               u.mile: '\\text{mi}', u.hbar: '\\hbar', u.coulomb_constant: '\\text{k}_e', u.angstrom: '\\text{Å}'} # Noen enheter har spesielle LaTeX-representasjoner som ikke kan hentes fra sympy, så de må legges inn manuelt her
    
    if enhet == ETTALLET:
        return ''
    if enhet in abbrevs:
        return abbrevs[enhet]
    if hasattr(enhet, '_latex_repr'):
        if enhet._latex_repr:
            return str(enhet._latex_repr)
    if '/' in str(enhet):
        teller, nevner = str(enhet).split('/')
        teller = latex_enhet(parse_enhet(teller)) if teller != '1' else '1'
        nevner = latex_enhet(parse_enhet(nevner)) if nevner != '1' else '1'
        if nevner != '1':
            return '\\frac{' + teller + '}{' + nevner + '}'
        else:
            return teller
    if enhet.is_Pow: # Sjekker om enheten er i form av x**(y)
        base = latex_enhet(enhet.base)
        eksponent = str(enhet.exp)
        if base != '':
            if '/' in eksponent:
                eksponent = '\\frac{' + eksponent.split('/')[0] + '}{' + eksponent.split('/')[1] + '}'
            return base + '^{' + eksponent + '}'
        else:
            return '1'
    if '*' in str(enhet):
        enhet = str(enhet).replace('**', '^')
        faktorer = str(enhet).split('*')
        latex_faktorer = [latex_enhet(parse_enhet(faktor)) for faktor in faktorer if (faktor != '' and faktor != '1')]
        return ' \\cdot '.join(latex_faktorer)
    
    if hasattr(enhet, 'abbrev'):
        return '\\text{' + str(enhet.abbrev) + '}'
    # Hvis ingen av de spesielle tilfellene gjelder, returneres enheten som den er, men i LaTeX-format
    return '\\text{' + str(enhet) + '}'
    

def sjekk_resultat(resultat: sp.Expr) -> tuple[sp.Rational, sp.Expr, sp.Mul]:
    """Denne funksjonen sjekker resultatet av en beregning og returnerer det rasjonale og irrasjonale delen og enheten til resultatet.
       Hvis resultatet er ugyldig, returneres det som NaN og en varseltekst blir satt.
    Args:
        resultat (sp.Expr): Tallet som skal sjekkes.

    Returns:
        sp.Rational:   Rasjonalt del av resultatet.
        sp.Expr:       Irrasjonalt del av resultatet.
        sp.Mul:        Enheten til resultatet.
    """
    global varseltekst

    forenklet_resultat = resultat.simplify()
    rasjonal, irrasjonal, målenhet = split_måltall_enhet(forenklet_resultat)

    if irrasjonal == IKKETALL:
        varseltekst = 'OBS! Ugyldig input. Resultatet er ikke et reelt tall.'
        rasjonal = IKKETALL
        irrasjonal = ETTALLET
        målenhet = ETTALLET
        return rasjonal, irrasjonal, målenhet
    if irrasjonal == NULLTALLET:
        rasjonal = NULLTALLET
        irrasjonal = ETTALLET
        målenhet = ETTALLET
        return rasjonal, irrasjonal, målenhet
    if irrasjonal.is_extended_negative:
        rasjonal *= -1
        irrasjonal *= -1
    if irrasjonal == sp.S.Infinity:
        rasjonal *= sp.S.Infinity
        irrasjonal = ETTALLET
    return rasjonal, irrasjonal, målenhet


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
        self.rasjonal: sp.Rational = ETTALLET   # Den rasjonale delen av tallet
        self.irrasjonal: sp.Expr = ETTALLET  # Den irrasjonale delen av tallet (Det er 1 hvis tallet er rasjonalt)
        self.feil: sp.Integer = NULLTALLET  # Feil i beregningen
        negativ: bool = False  # Er tallet negativt?
        self.enhet: sp.Mul = ETTALLET # Enhet for tallet, f.eks. u.meter for meter
        self.latexenhet: str = '' # LaTeX-representasjon av enheten, f.eks. 'm' for meter

        # Det fungerer dessverre ikke med match case. Derfor må jeg bruke if.
        if tall == '': # Det er ikke et tall
            return
        
        # Sjekker om tallet har en enhet
        if '  ' in tall:
            tall, enhet = tall.split('  ')
            self.enhet = parse_enhet(enhet)
            self.rasjonal, self.irrasjonal, self.enhet = split_måltall_enhet(self.enhet)
            self.latexenhet = latex_enhet(self.enhet)
        
        # Sjekker om tallet er negativt eller positivt, og fjerner fortegnet
        if tall[0] == '-':
            negativ = True
            tall = tall[1:]
        elif tall[0] == '+':
            tall = tall[1:]
        
        if tall == '0/0': # Det er ikke et tall
            return       
        if tall == '0':
            self.rasjonal = NULLTALLET
            return
        if tall == '1':
            self.rasjonal *= sp.S.NegativeOne if negativ else ETTALLET
            return
        if tall == '0,5' or tall == '1/2':
            self.rasjonal *= sp.S.Half
            if negativ:
                self.rasjonal *= sp.S.NegativeOne
            return
        if tall == UENDELIG or tall == '1/0' or tall == 'oo':
            self.rasjonal = sp.S.NegativeInfinity if negativ else sp.S.Infinity
            return
        if tall == PI or tall == 'pi':
            self.rasjonal *= sp.S.NegativeOne if negativ else ETTALLET
            self.irrasjonal *= sp.S.Pi
            return
        if tall == TAU or tall == 'tau': # τ = 2π
            self.rasjonal *= sp.Integer(-2) if negativ else sp.Integer(2)
            self.irrasjonal *= sp.S.Pi
            return
        if tall == 'e':
            self.rasjonal *= sp.S.NegativeOne if negativ else ETTALLET
            self.irrasjonal *= sp.S.Exp1
            return
        if tall == FI or tall == 'fi' or tall == 'phi': # Det gyldne snitt φ
            self.rasjonal *= sp.S.NegativeOne if negativ else ETTALLET
            self.irrasjonal *= sp.S.GoldenRatio
            return
        if tall == 'c': # Speed of light in vacuum
            self.rasjonal = sp.S.NegativeOne if negativ else ETTALLET
            self.enhet = u.c
            self.latexenhet = 'c'
            return
        if tall == 'G': # Gravitational constant
            self.rasjonal = sp.S.NegativeOne if negativ else ETTALLET
            self.enhet = u.G
            self.latexenhet = 'G'
            return
        if tall == 'ħ': # Planck constant
            self.rasjonal = sp.S.NegativeOne if negativ else ETTALLET
            self.enhet = u.hbar
            self.latexenhet = '\\hbar'
            return
        if tall == 'k': # Coulomb constant
            self.rasjonal = sp.S.NegativeOne if negativ else ETTALLET
            self.enhet = u.coulomb_constant
            self.latexenhet = 'k_e'
            return
        if tall == 'R': # Gas constant
            self.rasjonal = sp.S.NegativeOne if negativ else ETTALLET
            self.enhet = u.R
            self.latexenhet = 'R'
            return
        if re.fullmatch(r'^\d+$', tall): # heltall
            self.rasjonal *= sp.Integer(tall)
            if negativ:
                self.rasjonal *= -1
            return
        if re.fullmatch(r'^\d+ \d+/\d+$', tall): # blandet tall
            hel, brøkdel = tall.split(' ')
            teller, nevner = brøkdel.split('/')
            self.rasjonal *= sp.Rational(sp.Integer(hel) + sp.Rational(teller, nevner))
            if negativ:
                self.rasjonal *= sp.S.NegativeOne
            return
        if re.fullmatch(r'^\d+/\d+$', tall): # brøk
            teller, nevner = tall.split('/')
            self.rasjonal = sp.Rational(teller, nevner)
            if negativ:
                self.rasjonal *= sp.S.NegativeOne
            return
        if re.fullmatch(r'^\d+,\d+$', tall): # desimaltall
            self.rasjonal *= sp.Rational(tall.replace(',', '.'))
            if negativ:
                self.rasjonal *= sp.S.NegativeOne
            return
        if re.fullmatch(r'^\d+(,\d+)?e[+-]?\d+$', tall): # vitenskapelig format
            self.rasjonal *= sp.Rational(tall.replace(',', '.'))
            if negativ:
                self.rasjonal *= sp.S.NegativeOne
            return


    def kopi(self):
        """Lager en kopi av et tall og dets enhet

        Returns:
            tall: et nytt tall med samme verdi som self
        """        
        c = tall('1')
        c.rasjonal = self.rasjonal
        c.irrasjonal = self.irrasjonal
        c.enhet = self.enhet
        c.latexenhet = self.latexenhet
        return c


    def gjør_hel(self) -> None:
        """Tar bare heltallig delen av et tall, dvs. 2,5 blir 2 og -2,5 blir -2.
           Enheten til tallet blir også beholdt.
        """
        if self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity: # ∞ eller -∞
            return
        if self.rasjonal != IKKETALL:
            self.rasjonal = sp.Integer(sp.Mul(self.rasjonal, self.irrasjonal))
            self.irrasjonal = ETTALLET 


    def resiprok(self) -> None:
        """Beregner resiprokverdien av et tall, dvs. 1/self
        """        
        if self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity: # 1/∞ = 0
            self.rasjonal = NULLTALLET
            self.irrasjonal = ETTALLET
            self.enhet = ETTALLET
            self.latexenhet = ''
            return
        
        if self.enhet != ETTALLET:
            self.enhet = ETTALLET / self.enhet
            self.latexenhet = latex_enhet(self.enhet)

        if self.rasjonal == NULLTALLET: # 1/0 = ∞ 
            self.rasjonal = sp.S.Infinity
            self.irrasjonal = ETTALLET
            self.enhet = ETTALLET
            self.latexenhet = ''
            return
        if self.rasjonal == ETTALLET and self.irrasjonal == ETTALLET: # 1/1 = 1
            return
        if self.rasjonal == sp.S.NegativeOne and self.irrasjonal == ETTALLET: # 1/(-1) = -1
            return
        if self.rasjonal != IKKETALL:
            b = self.kopi()
            self.rasjonal = ETTALLET/b.rasjonal
            self.irrasjonal = ETTALLET/b.irrasjonal
            return


    def pluss(self, addend) -> None:
        """Plusser sammen to tall: self + addend

        Args:
            addend (tall): Tallet som legges til
        """
        global varseltekst

        if addend.__class__ != tall: # Kan bare addere et tall
            return
        if addend.rasjonal == IKKETALL: # kan bare addere tall
            return
        if self.rasjonal == IKKETALL: # NaN + x = x
            self.rasjonal = addend.rasjonal
            self.irrasjonal = addend.irrasjonal
            self.enhet = addend.enhet
            self.latexenhet = addend.latexenhet
            return
        if self.rasjonal == sp.S.NegativeInfinity and addend.rasjonal == sp.S.Infinity: # -∞ + ∞ er ikke definert
            varseltekst = '-' + UENDELIG + ' pluss ' + UENDELIG + ' er ikke definert.'
            self.rasjonal = IKKETALL
            self.irrasjonal = ETTALLET
            self.enhet = ETTALLET
            self.latexenhet = ''
            return
        if self.rasjonal == sp.S.Infinity and addend.rasjonal == sp.S.NegativeInfinity: # ∞ - ∞ er ikke definert
            varseltekst = UENDELIG + ' minus ' + UENDELIG + ' er ikke definert.'
            self.rasjonal = IKKETALL
            self.irrasjonal = ETTALLET
            self.enhet = ETTALLET
            self.latexenhet = ''
            return
        if self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity: # ∞ + x = ∞
            return
        if addend.rasjonal == sp.S.Infinity or addend.rasjonal == sp.S.NegativeInfinity: # x + ∞ = ∞
            self.rasjonal = addend.rasjonal
            self.irrasjonal = addend.irrasjonal
            self.enhet = addend.enhet
            self.latexenhet = addend.latexenhet
            return
        
        if self.enhet != ETTALLET or addend.enhet != ETTALLET: # Hvis en av tallene har en enhet, må de ha samme enhet for å kunne adderes
            if not samme_dimensjon(self.enhet, addend.enhet):
                varseltekst = 'OBS! Kan ikke addere tall med forskjellige dimensjoner.'
                return
            if self.enhet != addend.enhet: # Hvis enhetene er forskjellige, må den ene konverteres til den andre før de kan adderes
                ny_enhet = u.convert_to(self.enhet, addend.enhet)
                rasjonal_måltall, irrasjonal_måltall, ny_enhet = split_måltall_enhet(ny_enhet)
                self.rasjonal *= rasjonal_måltall
                self.irrasjonal *= irrasjonal_måltall
                self.enhet = ny_enhet
        
        if self.irrasjonal == ETTALLET and addend.irrasjonal == ETTALLET: # Det er to rasjonale tall
            self.rasjonal += addend.rasjonal
            self.latexenhet = addend.latexenhet
            return

        verdi = sp.Add(sp.Mul(self.rasjonal, self.irrasjonal, self.enhet), sp.Mul(addend.rasjonal, addend.irrasjonal, addend.enhet))
        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(verdi)
        self.latexenhet = latex_enhet(self.enhet)

    
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
        if self.rasjonal == IKKETALL: # kan bare gange tall
            varseltekst = 'Du kan bare gange tall.'
            return
        if faktor.rasjonal == IKKETALL: # kan bare gange med tall
            varseltekst = 'Du kan bare gange tall.'
            return
        
        if (self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity) \
            and (faktor.rasjonal == sp.S.Infinity or faktor.rasjonal == sp.S.NegativeInfinity): # ∞ * ∞ er ikke definert
            varseltekst = UENDELIG + ' ganger ' + UENDELIG + ' er ikke definert.'
            self.rasjonal = IKKETALL
            self.enhet = ETTALLET
            self.latexenhet = ''
            return
        if self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity: # ∞ * x er ∞
            match faktor.rasjonal:
                case sp.S.Zero:
                    varseltekst = UENDELIG + ' ganger 0 er ikke definert.'
                    self.rasjonal = IKKETALL
                    self.enhet = ETTALLET
                    self.latexenhet = ''
                case sp.S.Infinity:
                    varseltekst = UENDELIG + ' ganger ' + UENDELIG + ' er ikke definert.'
                    self.rasjonal = IKKETALL
                    self.enhet = ETTALLET
                    self.latexenhet = ''
                case sp.S.NegativeInfinity:
                    varseltekst = UENDELIG + ' ganger -' + UENDELIG + ' er ikke definert.'
                    self.rasjonal = IKKETALL
                    self.enhet = ETTALLET
                    self.latexenhet = ''
                case _:
                    self.rasjonal *= faktor.rasjonal
                    self.enhet *= faktor.enhet
                    self.latexenhet = latex_enhet(self.enhet)
            self.irrasjonal = ETTALLET
            return
        if faktor.rasjonal == sp.S.Infinity or faktor.rasjonal == sp.S.NegativeInfinity: # x * ∞ er ∞
            if self.rasjonal == NULLTALLET:
                varseltekst = '0 ganger ' + UENDELIG + ' er ikke definert.'
                self.rasjonal = IKKETALL
                self.enhet = ETTALLET
                self.latexenhet = ''
            else:
                self.rasjonal *= faktor.rasjonal
                self.enhet *= faktor.enhet
                self.latexenhet = latex_enhet(self.enhet)
            self.irrasjonal = ETTALLET
            return
        if self.rasjonal == NULLTALLET or faktor.rasjonal == NULLTALLET: # x * 0 er 0
            self.rasjonal = NULLTALLET
            self.irrasjonal = ETTALLET
            self.enhet = ETTALLET
            self.latexenhet = ''
            return

        resultat = sp.Mul(self.rasjonal, self.irrasjonal, self.enhet, faktor.rasjonal, faktor.irrasjonal, faktor.enhet)
        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)
        
        
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

        if self.irrasjonal != ETTALLET or self.rasjonal < NULLTALLET or not self.rasjonal.is_integer:
            varseltekst = 'OBS! Kan bare beregne fakultet av naturlige tall.'
            return
        self.rasjonal = sp.Integer(sp.factorial(self.rasjonal))

  
    def fib(self) -> None:
        """Beregner fibonacci-tallet av et naturlig tall
        """        
        global varseltekst

        if self.irrasjonal != ETTALLET or self.rasjonal < NULLTALLET or not self.rasjonal.is_integer:
            varseltekst = 'OBS! Kan bare beregne fibinacci-tallet av naturlige tall.'
            return
        self.rasjonal = sp.fibonacci(self.rasjonal)

  
    def binom(self, k) -> None:
        """Beregner binominalkoeffisienten self over k

        Args:
            k (tall): nederste tall i binominalkoeffisienten n over k
        """        
        global varseltekst

        if k.__class__ != tall: # Kan bare regne med tall
            return
        if k.irrasjonal != ETTALLET or not k.rasjonal.is_integer:
            varseltekst = 'OBS! Kan bare beregne binominalkoeffisienten over et heltall.'
            return
        if k.rasjonal < NULLTALLET:
            self.rasjonal = NULLTALLET
            self.irrasjonal = ETTALLET
            return
        if self.irrasjonal == ETTALLET:
            self.rasjonal = sp.binomial(self.rasjonal, k.rasjonal)
        else:
            varseltekst = 'OBS! Kan bare beregne binominalkoeffisienten av rasjonale tall.'

  
    def opphøyd_i(self, potens) -> None:
        """Beregner self opphøyd i potens, dvs. opphøyer et tall i en potens

        Args:
            potens (tall): potensen som tallet opphøyes i
        """        
        global varseltekst

        if potens.__class__ != tall or potens.enhet != ETTALLET: # Kan bare regne med tall
            return
        if potens.rasjonal == NULLTALLET:
            self.rasjonal = ETTALLET
            self.irrasjonal = ETTALLET
            self.enhet = ETTALLET
            self.latexenhet = ''
            return
        if potens.rasjonal == ETTALLET and potens.irrasjonal == ETTALLET:
            return

        resultat = sp.Pow(sp.Mul(self.rasjonal, self.irrasjonal, self.enhet), sp.Mul(potens.rasjonal, potens.irrasjonal))
        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)


    def lb(self) -> None:
        """Beregner logaritme til base 2
        """        
        global varseltekst

        if self.rasjonal < NULLTALLET:
            varseltekst = 'OBS! Kan bare beregne logaritmus av positive tall!'
            return
        if self.rasjonal == IKKETALL or self.enhet != ETTALLET:
            varseltekst = 'OBS! Kan bare beregne logaritmus av et tall!'
            return
        if self.rasjonal == NULLTALLET: # lb(0) = -∞
            self.rasjonal = sp.S.NegativeInfinity
            self.irrasjonal = ETTALLET
            return
        if self.rasjonal == sp.S.Infinity: # lb(∞) = ∞
            return
        
        # Beregn logaritmen til base 2 ved hjelp av sympy
        resultat = sp.log(sp.Mul(self.rasjonal, self.irrasjonal), 2)
        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)
        

    def lg(self) -> None:
        """Beregner logaritme til base 10
        """        
        global varseltekst

        if self.rasjonal < NULLTALLET:
            varseltekst = 'OBS! Kan bare beregne logaritmus av positive tall!'
            return
        if self.rasjonal == IKKETALL or self.enhet != ETTALLET:
            varseltekst = 'OBS! Kan bare beregne logaritmus av et tall!'
            return
        if self.rasjonal == NULLTALLET: # lg(0) = -∞
            self.rasjonal = sp.S.NegativeInfinity
            self.irrasjonal = ETTALLET
            return
        if self.rasjonal == sp.S.Infinity: # lg(∞) = ∞
            return
        
        # Beregn logaritmen til base 10 ved hjelp av sympy
        resultat = sp.log(sp.Mul(self.rasjonal, self.irrasjonal), 10)
        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)

    def ln(self) -> None:
        """Beregner logaritme til base e, den naturlige logaritme
        """        
        global varseltekst

        if self.rasjonal < NULLTALLET:
            varseltekst = 'OBS! Kan bare beregne logaritmus av positive tall!'
            return
        if self.rasjonal == IKKETALL or self.enhet != ETTALLET:
            varseltekst = 'OBS! Kan bare beregne logaritmus av et tall!'
            return
        if self.rasjonal == NULLTALLET: # ln(0) = -∞
            self.rasjonal = sp.S.NegativeInfinity
            self.irrasjonal = ETTALLET
            return
        if self.rasjonal == sp.S.Infinity: # ln(∞) = ∞
            return
        
        # Beregn logaritmen til base e ved hjelp av sympy
        resultat = sp.log(sp.Mul(self.rasjonal, self.irrasjonal), sp.S.Exp1)
        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)


    def sin(self) -> None:
        """Beregner sinus av en vinkel i grader eller radianer
        """        
        global varseltekst

        if self.rasjonal == IKKETALL:
            varseltekst = 'OBS! Kan bare beregne sinus av et tall!'
            return
        if self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity:
            varseltekst = 'OBS! Kan ikke beregne sinus av uendelig!'
            self.rasjonal = IKKETALL
            return
        
        # Beregn sinus ved hjelp av sympy
        if self.enhet == u.degree or self.enhet == u.steradian:
            ny_enhet = u.convert_to(self.enhet, u.radian)
            rasjonal, irrasjonal, ny_enhet = split_måltall_enhet(ny_enhet)
            resultat = sp.sin(sp.Mul(self.rasjonal, self.irrasjonal, rasjonal, irrasjonal))
        elif self.enhet == u.radian or self.enhet == ETTALLET:
            resultat = sp.sin(sp.Mul(self.rasjonal, self.irrasjonal))
        else:
            varseltekst = 'OBS! Kan bare beregne sinus av vinkler i grader eller radianer!'
            self.rasjonal = IKKETALL
            return

        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)


    def cos(self) -> None:
        """Beregner kosinus av en vinkel i grader eller radianer
        """        
        global varseltekst

        if self.rasjonal == IKKETALL:
            varseltekst = 'OBS! Kan bare beregne kosinus av et tall!'
            return
        if self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity:
            varseltekst = 'OBS! Kan ikke beregne kosinus av uendelig!'
            self.rasjonal = IKKETALL
            return
        
        # Beregn kosinus ved hjelp av sympy
        if self.enhet == u.degree or self.enhet == u.steradian:
            ny_enhet = u.convert_to(self.enhet, u.radian)
            rasjonal, irrasjonal, ny_enhet = split_måltall_enhet(ny_enhet)
            resultat = sp.cos(sp.Mul(self.rasjonal, self.irrasjonal, rasjonal, irrasjonal))
        elif self.enhet == u.radian or self.enhet == ETTALLET:
            resultat = sp.cos(sp.Mul(self.rasjonal, self.irrasjonal))
        else:
            varseltekst = 'OBS! Kan bare beregne kosinus av vinkler i grader eller radianer!'
            self.rasjonal = IKKETALL
            return
        
        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)


    def tan(self) -> None:
        """Beregner tangens av en vinkel i grader eller radianer
        """        
        global varseltekst

        if self.rasjonal == IKKETALL:
            varseltekst = 'OBS! Kan bare beregne tangens av et tall!'
            return
        if self.rasjonal == sp.S.Infinity or self.rasjonal == sp.S.NegativeInfinity:
            varseltekst = 'OBS! Kan ikke beregne tangens av uendelig!'
            self.rasjonal = IKKETALL
            return
        
        # Beregn tangens ved hjelp av sympy
        if self.enhet == u.radian or self.enhet == ETTALLET:
            if self.irrasjonal == sp.S.Pi:
                if abs(self.rasjonal) % ETTALLET == sp.S.Half:
                    varseltekst = 'OBS! Tangens av π/2 + nπ er ikke definert!'
                    self.rasjonal = IKKETALL
                    self.irrasjonal = ETTALLET
                    return
            resultat = sp.tan(sp.Mul(self.rasjonal, self.irrasjonal))
        elif self.enhet == u.degree or self.enhet == u.steradian:
            ny_enhet = u.convert_to(self.enhet, u.radian)
            rasjonal, irrasjonal, ny_enhet = split_måltall_enhet(ny_enhet)
            argument = sp.Mul(self.rasjonal, self.irrasjonal, rasjonal, irrasjonal)
            if abs(argument) % sp.Integer(180) == 90:
                varseltekst = 'OBS! Tangens av 90° + n*180° er ikke definert!'
                self.rasjonal = IKKETALL
                self.irrasjonal = ETTALLET
                return
            resultat = sp.tan(argument)
        else:
            varseltekst = 'OBS! Kan bare beregne tangens av vinkler i grader eller radianer!'
            self.rasjonal = IKKETALL
            return

        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.latexenhet = latex_enhet(self.enhet)


    def arcsin(self) -> None:
        """Beregner arcussinus av et tall i intervall [-1, 1]
        """        
        global varseltekst

        if self.rasjonal == IKKETALL or self.enhet != ETTALLET:
            varseltekst = 'OBS! Kan bare beregne arcsinus av et tall!'
            return
        if abs(sp.Mul(self.rasjonal, self.irrasjonal)) > ETTALLET:
            varseltekst = 'OBS! Tallet må være i intervall [-1, 1]!'
            self.rasjonal = IKKETALL
            return
        
        # Beregn arcsinus ved hjelp av sympy
        resultat = sp.asin(sp.Mul(self.rasjonal, self.irrasjonal))

        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.omgjør(u.degree) # Omgjør resultatet til grader


    def arccos(self) -> None:
        """Beregner arcuscosinus av et tall i intervall [-1, 1]
        """        
        global varseltekst
        if self.rasjonal == IKKETALL or self.enhet != ETTALLET:
            varseltekst = 'OBS! Kan bare beregne arccosinus av et tall!'
            return
        if abs(sp.Mul(self.rasjonal, self.irrasjonal)) > ETTALLET:
            varseltekst = 'OBS! Tallet må være i intervall [-1, 1]!'
            self.rasjonal = IKKETALL
            return
        
        # Beregn arccosinus ved hjelp av sympy
        resultat = sp.acos(sp.Mul(self.rasjonal, self.irrasjonal))

        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.omgjør(u.degree) # Omgjør resultatet til grader


    def arctan(self) -> None:
        """Beregner arctangens av et tall
        """        
        global varseltekst
        if self.rasjonal == IKKETALL or self.enhet != ETTALLET:
            varseltekst = 'OBS! Kan bare beregne tangens av et tall!'
            return
        
        # Beregn arctangens ved hjelp av sympy
        resultat = sp.atan(sp.Mul(self.rasjonal, self.irrasjonal))

        self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(resultat)
        self.omgjør(u.degree) # Omgjør resultatet til grader


    def utvid(self) -> None:
        """Utvider et uttrykk ved å bruke sympys expand(func=True)
        """
        if self.rasjonal != IKKETALL and self.irrasjonal != ETTALLET:
            uttrykk = sp.expand(sp.Mul(self.rasjonal, self.irrasjonal, self.enhet), func=True)
            self.rasjonal, self.irrasjonal, self.enhet = sjekk_resultat(uttrykk)
            self.latexenhet = latex_enhet(self.enhet)


    def evaluer(self) -> None:
        """Estimerer et irrasjonalt tall til et rasjonalt tall ved å bruke sympy
        """
        if self.rasjonal != IKKETALL and self.irrasjonal != ETTALLET:
            tall = sp.Mul(self.rasjonal, self.irrasjonal)           # Gange den rasjonale med den irrasjonale delen
            rasjonal = sp.Rational(str(tall.evalf(n=MAKS_SIFRE)))   # Evaluer og omgjør til et rasjonalt tall
            self.rasjonal, _, _ = sjekk_resultat(rasjonal)          # Sjekk resultatet
            self.irrasjonal = ETTALLET
            avvik = abs(sp.Add(tall, -self.rasjonal)).evalf()       # Beregne avviket for å beregne feilen
            self.feil = NULLTALLET if avvik < 1e-100 else sp.Integer(sp.log(avvik, 10).evalf())


    def omgjør(self, ny_enhet) -> None:
        """Omgjør et tall til en annen enhet ved å gange tallet med konverteringsfaktoren mellom de to enhetene

        Args:
            ny_enhet: enheten som tallet skal omgjøres til
        """
        global varseltekst

        if self.rasjonal == IKKETALL:
            varseltekst = 'OBS! Kan bare omgjøre tall.'
            return
        if self.enhet == ETTALLET:
            # Hvis tallet inneholder π antas at det er radianer som skal omgjøres til grader
            if ny_enhet == u.degree and self.irrasjonal == sp.S.Pi:
                self.rasjonal *= 180
                self.irrasjonal = ETTALLET
            self.enhet = ny_enhet
            self.latexenhet = latex_enhet(self.enhet)
            return
        if ny_enhet == ETTALLET: # Det vil aldri skje siden omgjør(ny_enhet) kalles kun med en ny enhet
            varseltekst = 'OBS! Kan ikke omgjøre til enheten 1.'
            return
        if not samme_dimensjon(self.enhet, ny_enhet):
            varseltekst = 'OBS! Kan ikke omgjøre til enheter med forskjellige dimensjoner.'
            return
        
        konverteringsfaktor = u.convert_to(self.enhet, ny_enhet)
        rasjonal_konverteringsfaktor, irrasjonal_konverteringsfaktor, _ = split_måltall_enhet(konverteringsfaktor)
        self.rasjonal *= rasjonal_konverteringsfaktor
        self.irrasjonal *= irrasjonal_konverteringsfaktor
        self.enhet = ny_enhet
        self.latexenhet = latex_enhet(self.enhet)


    def omgjør_SI(self) -> None:
        """Omgjør et tall til SI-enhet ved å gange tallet med konverteringsfaktoren mellom den opprinnelige enheten og SI-enheten
        """
        global varseltekst

        if self.rasjonal == IKKETALL:
            varseltekst = 'OBS! Kan bare omgjøre tall.'
            return
        if self.enhet == ETTALLET: # Hvis tallet ikke har en enhet, kan det ikke omgjøres
            return
        
        konverteringsfaktor = u.convert_to(self.enhet, u.systems.SI._base_units)
        rasjonal_konverteringsfaktor, irrasjonal_konverteringsfaktor, enhet = split_måltall_enhet(konverteringsfaktor)
        self.rasjonal *= rasjonal_konverteringsfaktor
        self.irrasjonal *= irrasjonal_konverteringsfaktor
        self.enhet = enhet
        self.latexenhet = latex_enhet(self.enhet)


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
        case 'rydd':
            tallene.clear()
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
                tallene[-1].rasjonal *= sp.S.NegativeOne # Bytter fortegn
        case '--': # bytt fortegn
            tallene[-1].rasjonal *= sp.S.NegativeOne # Bytter fortegn
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
        case 'log': # logaritme til base 10
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
        case 'sum':
            while len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].pluss(a)
        case 'Ø':
            antall = tall(str(len(tallene)))
            while len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].pluss(a)
            tallene[-1].delt_med(antall)
        case 'prod':
            while len(tallene) > 1:
                a = tallene.pop()
                tallene[-1].ganger(a)
        case 'hel':
            tallene[-1].gjør_hel()
        case 'utvid': # utvider et uttrykk
            if tallene[-1].irrasjonal != ETTALLET:
                tallene[-1].utvid()
        case 'eval': # evaluerer et uttrykk
            if tallene[-1].irrasjonal != ETTALLET:
                tallene[-1].evaluer()
        case 'SI': # omgjør til SI-enhet
            tallene[-1].omgjør_SI()
        case _:
            enhet = parse_enhet(operasjon)
            if enhet != ETTALLET:
                tallene[-1].omgjør(enhet)
            else:
                varseltekst = 'OBS! Ukjent operasjon'


def notasjon(tallet: sp.Rational, feil: sp.Integer) -> tuple[str, sp.Integer]:
    """Konverterer en brøk til vitenskapelig notasjon

    Args:
        tallet (sp.Rational): tallet som skal skrives i vitenskapelig notasjon
        feil (sp.Integer): feil i beregningen
    Returns:
        str: tallet i vitenskapelig notasjon slik at den kan brukes i LaTeX
        sp.Integer: den oppdaterte feilen
    """
    if not isinstance(tallet, sp.Rational):
        return '', feil
    fortegn: str = '-' if tallet < NULLTALLET else ''
    ellipsis: str = ''
    teller: sp.Integer = abs(sp.numer(tallet))
    nevner: sp.Integer = sp.denom(tallet)
    scientific_str: str = ''

    if teller < nevner:
        # Beregn hvor mye mindre enn 1 tallet er, dvs. hvor mange nuller bak kommaet det har
        rest: sp.Integer = teller
        nuller: sp.Integer = NULLTALLET

        # Tell hvor ofte vi må gange rest med ti til det blir større enn nevner
        while rest < nevner:
            rest *= sp.Integer(10)
            nuller += ETTALLET
        
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
            rest *= sp.Integer(10) 
            indeks += ETTALLET
        
        if len(desimaldel) >= MAKS_SIFRE: # Vi har nådd ønsket presisjon, men ikke funnet en periode
            førperiode = ''.join(desimaldel)
            periode = ''
            feil = max(-MAKS_SIFRE, feil) if feil != NULLTALLET else -MAKS_SIFRE # Feilen er minst 10^-MAKS_SIFRE
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
                nuller -= ETTALLET # alt flyttes en plass til venstre
            periode = ''

        if nuller == 1: # Det er bare en null foran kommaet
            if periode:
                scientific_str = '0,' + førperiode + '\\overline{' + periode + '}'
            elif førperiode:
                scientific_str = '0,' + førperiode + ellipsis
            else:
                scientific_str = '0'
        elif nuller < MAKS_SIFRE//sp.Integer(2): # Hvis det er ikke altfor mange nuller, kan tallet framstilles som desimaltall
            if periode:
                if førperiode:
                    scientific_str = '0,' + ('0' * int(nuller - 1)) + førperiode + '\\overline{' + periode + '}'
                    # OBS: Det er (nuller - 1) fordi en null står foran kommaet
                elif periode[-1] == '0': # Flytt nuller fra slutten av perioden til begynnelsen
                    uten_nuller = periode.rstrip('0') # Det kan være flere 0er
                    n_nuller = len(periode) - len(uten_nuller) # Beregn hvor mange nuller det er
                    scientific_str = '0,' + ('0' * int(nuller - n_nuller - ETTALLET)) + '\\overline{' + ('0' * n_nuller) + uten_nuller + '}'
                else: # Det er ikke nuller som må flyttes
                    scientific_str = '0,' + ('0' * int(nuller - 1)) + '\\overline{' + periode + '}'
            elif førperiode:
                scientific_str = '0,' + ('0' * int(nuller - 1)) + førperiode + ellipsis
            else:
                scientific_str = '0'
        else: # Tallet må framstilles i vitenskapelig notasjon
            if feil != NULLTALLET:
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
        heldel = sp.Integer(abs(tallet))
        helstr = str(heldel)
        eksponent = len(helstr) - 1  # eksponenten er antall sifre i heldelen minus 1 (minus 1 fordi det skal være ett siffer foran kommaet)
        # Beregn desimaldelen ved lang divisjon
        rest = teller % nevner * sp.Integer(10)
        huskelapp: dict = {}  # rest -> indeks i desimaldelen
        desimaldel: list[str] = []
        indeks: int = 0
        # Beregn sifrene til vi finner perioden eller når ønsket presisjon
        while (rest not in huskelapp) and (len(desimaldel) < MAKS_SIFRE):
            siffer = rest // nevner
            desimaldel.append(str(siffer))
            huskelapp[rest] = indeks # Husk at vi fant rest ved indeks
            rest %= nevner # Beregne ny rest
            rest *= sp.Integer(10)
            indeks += ETTALLET
            
        if len(desimaldel) >= MAKS_SIFRE: # Vi har nådd ønsket presisjon, men ikke funnet en periode
            førperiode = ''.join(desimaldel)
            periode = ''
            feil = max(-MAKS_SIFRE, feil) if feil != NULLTALLET else -MAKS_SIFRE # Feilen er minst 10^-MAKS_SIFRE
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
                eksponent -= ETTALLET # alt flyttes en plass til venstre
            periode = ''

        if eksponent <= ETTALLET: # Hvis heldelen har bare en eller to sifre, kan de bli stående
            første_siffer = helstr
            brøkdel = førperiode + '\\overline{' + periode + '}' if periode else førperiode
            if brøkdel:
                scientific_str = første_siffer + ',' + brøkdel
            else:
                scientific_str = første_siffer
        elif eksponent < MAKS_SIFRE//sp.Integer(4): # Hvis tallet ikke er for stor, kan den framstilles som desimaltall
            første_siffer = '{:,}'.format(int(heldel)).replace(',', '~')
            brøkdel = førperiode + '\\overline{' + periode + '}' if periode else førperiode
            eksponent = NULLTALLET
            if brøkdel:
                scientific_str: str = første_siffer + ',' + brøkdel
            else:
                scientific_str: str = første_siffer
        else: # ellers bruker vi vitenskapelig notasjon med ett siffer foran kommaet
            if feil != NULLTALLET:
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

    if feil != NULLTALLET:
        return '\\approx ' + fortegn + scientific_str, feil
    return '= ' + fortegn + scientific_str, feil 


def heltall_til_latex(verdi: sp.Integer, feil: sp.Integer) -> str:
    """Konverterer et heltall til en LaTeX-streng.
    Args:
        verdi (sp.Integer): Heltallet som skal konverteres til LaTeX-format.
        feil (sp.Integer): Feil i beregningen
    Returns:
        str: LaTeX-strengen som representerer tallet
    """
    likhetstegn: str = '= '
    ellipsis: str = ''

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
            feil = max(eksponent - MAKS_SIFRE - ikke_bare_nuller - 1, feil) if feil != NULLTALLET else eksponent - MAKS_SIFRE - ikke_bare_nuller - 1
        resten_av_tall = ''.join(str(sifr) for sifr in sifre[2:]) 

    uten_nuller = resten_av_tall.rstrip('0') # Fjerner nuller på slutten av resten av tallet

    komma = ',' if len(uten_nuller) > 0 else '' # Hvis det er ingen sifre etter kommaet, fjerner vi kommaet

    return likhetstegn + første_siffer + komma + uten_nuller + ellipsis + '\\cdot 10^{' + str(eksponent) + '}', feil


def til_latex(tallet: tall) -> str:
    """Konverterer et tall til en LaTeX-streng.
    Args:
        tallet (tall): Tallet som skal konverteres til LaTeX-format

    Returns:
        str: LaTeX-strengen som representerer tallet
    """
    if tallet.rasjonal == sp.S.Infinity:
        return '= \\infty'
    if tallet.rasjonal == sp.S.NegativeInfinity:
        return '= -\\infty'
    if tallet.rasjonal == IKKETALL:
        return '\\bot'

    verdi = tallet.rasjonal if isinstance(tallet.rasjonal, sp.Rational) else None

    limit = sp.Pow(10, MAKS_SIFRE)
    if verdi and tallet.irrasjonal == ETTALLET: # Rasjonale tall kan vises som som brøk, blandet tall eller heltall
        likhetstegn: str = '= ' if tallet.feil == NULLTALLET else '\\approx '
        if verdi.is_integer: # Det er et heltall
            if abs(verdi) < limit:
                # Formater tallet slik at de er lettere å lese, dvs. 1000000 blir 1 000 000
                return likhetstegn + '{:,}'.format(int(verdi)).replace(',', '~')
            # Et heltall med flere enn MAKS_SIFRE sifre vises i vitenskapelig notasjon
            uttrykk, feil = heltall_til_latex(sp.Integer(verdi), tallet.feil)
            tallet.feil = feil
            return uttrykk
        if verdi.is_rational: # Det er en brøk
            teller: sp.Integer = sp.numer(verdi) # Henter teller fra brøken
            nevner: sp.Integer = sp.denom(verdi) # Henter nevner fra brøken
            if verdi > 1 or verdi < -1:
                if max(teller, nevner) < limit: # Det er et blandet tall som kan vises på vanlig måte og vitenskapelig notasjon
                    hel = sp.Integer(verdi)
                    # Formater hel, teller og nevner slik at de er lettere å lese, dvs. 1000000 blir 1 000 000
                    nevnerstr = '{:,}'.format(int(nevner)).replace(',', '~')
                    tellerstr = '{:,}'.format(int(abs(sp.Add(teller, -sp.Mul(hel, nevner))))).replace(',', '~')
                    helstr = likhetstegn + '{:,}'.format(int(hel)).replace(',', '~')
                    uttrykk, feil = notasjon(verdi, tallet.feil)
                    tallet.feil = feil
                    return helstr + ' \\frac{' + tellerstr + '}{' + nevnerstr + '} ' + uttrykk
                uttrykk, feil = notasjon(verdi, tallet.feil)
                tallet.feil = feil
                return uttrykk  
            teller = abs(teller)
            if max(teller, nevner) < limit: # Det er en brøk som kan vises på vanlig måte
                # Formater teller og nevner slik at de er lettere å lese, dvs. 1000000 blir 1 000 000
                tellerstr = '{:,}'.format(int(teller)).replace(',', '~')
                nevnerstr = '{:,}'.format(int(nevner)).replace(',', '~')
                if verdi < NULLTALLET: # Hvis tallet er negativt, må vi vise det med minus foran
                    uttrykk, feil = notasjon(-verdi, tallet.feil)
                    tallet.feil = feil
                    return likhetstegn + '-\\frac{' + tellerstr + '}{' + nevnerstr + '}' + uttrykk
                else:
                    uttrykk, feil = notasjon(verdi, tallet.feil)
                    tallet.feil = feil
                    return likhetstegn + '\\frac{' + tellerstr + '}{' + nevnerstr + '}' + uttrykk
            # Hvis teller eller nevner er for store, brukes bare vitenskapelig notasjon
            uttrykk, feil = notasjon(verdi, tallet.feil)
            tallet.feil = feil
            return uttrykk
        return '\\bot' # Det er ikke et tall som kan formateres til LaTeX-formatet

    # Tallet er irrasjonalt
    return '= ' + sp.latex(sp.Mul(verdi, tallet.irrasjonal), decimal_separator='comma', max=MAKS_SIFRE)


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
    uttrykk_tall = '$' + til_latex(tallet) + '$'
    uttrykk_feil = '$(\\pm10^{' + str(tallet.feil) + '})$' if tallet.feil != NULLTALLET else ''
    uttrykk_enhet = '$' + tallet.latexenhet + '$' if tallet.latexenhet else ''
       
    # Tilføy uttrykket som tekst til figuren
    if uttrykk_feil:
        ax.text(0.45, 0.4, uttrykk_tall, transform=fig.transFigure, color='black', fontsize=16, ha='center', va='baseline')
        ax.text(0.90, 0.4, uttrykk_feil, transform=fig.transFigure, color='gray', fontsize=14, ha='center', va='baseline')
    else:
        ax.text(0.5, 0.4, uttrykk_tall, transform=fig.transFigure, color='black', fontsize=16, ha='center', va='baseline')
    if uttrykk_enhet:
        ax.text(0.96, 0.4, uttrykk_enhet, transform=fig.transFigure, color='black', fontsize=16, ha='center', va='baseline')

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
    if len(stabel) > 0:
        resultat = stabel[-1].kopi()
        resultat.evaluer()
        img = skriv_resultat(resultat)
        lbl_resultat.image = toga.Image(img)


class SymKalkulator(toga.App):
    """Det er klassen som lager app'en 'SymKalkulator'.

    Args:
        toga (App): Toga er GUI'en som brukes
    """    
    def startup(self) -> None:
        """Startup er en funksjon som hver Toga-app må ha. Det er funksjonen som kjøres når SymKalkulator() kalles.
        """        
        global varseltekst, spesielle_tall
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
            elif innput in spesielle_tall: # Det har blitt tastet inn flere enn ett tegn
                nytall = tall(innput)
                stabel.append(nytall)
            elif len(innput) == 1: # Ikke ett tall, men bare ett tegn har blitt tastet inn
                # Det er en operasjon fra listen ['+', '-', '*', '/', '%', '!', 'v']
                operasjon = innput
                beregne(stabel, operasjon) # Anvender operasjonen på tallene som ligger på stabelen
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
            elif '  ' in innput: # Det er et tall med en enhet, f.eks. '5 m'
                nytall = tall(innput)
                stabel.append(nytall)
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
                pyperclip.copy(str(sp.Mul(tallet.rasjonal, tallet.irrasjonal).evalf(n=MAKS_SIFRE)))


        # Det følgende er fortsatt del av startup-funksjonen. Her defineres hvordan vinduet ser ut og fungerer.
        innhold = toga.Box(style=Pack(direction=COLUMN, margin=2)) # Innholdet i hovedvinduet
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
        enheter_cmd = toga.Command(
            self.enheter,
            text='Enheter',
            tooltip='Vis enheter du kan bruke',
            shortcut=toga.Key.MOD_1 + 'e',
            group=toga.Group.HELP
        )
        copy_cmd = toga.Command(
            kopier,
            text='Kopier',
            tooltip='Kopier resultat til clipboard',
            shortcut=toga.Key.MOD_1 + 'c'
        )
        self.commands.add(hjelp_cmd)
        self.commands.add(enheter_cmd)
        self.commands.add(copy_cmd)
        
        # Lage stabel-boxen
        stabel_box = toga.Box(style=Pack(direction=COLUMN, align_items=CENTER, margin=2))
        global lbl_stabel
        lbl_stabel = [toga.ImageView(None, style=Pack(height=80, margin=2)) for _ in range(8)]
        for i in range(7):
            stabel_box.add(lbl_stabel[i])
        stabel_box.add(toga.Divider())
        stabel_box.add(lbl_stabel[7])

        # Lage innput-boxen
        input_box = toga.Box(style=Pack(direction=COLUMN, align_items=CENTER, margin=2))
        self.innput = toga.TextInput(style=Pack(text_align=CENTER, width=1200, height=50, margin=2, font_size=20), on_confirm=inntastet)
        input_box.add(self.innput)

        # Lage resultat-boxen
        resultat_box = toga.Box(style=Pack(direction=COLUMN, align_items=CENTER, margin=2))
        global lbl_resultat
        lbl_resultat = toga.ImageView(None, style=Pack(height=80, margin=2))
        resultat_box.add(lbl_resultat)

        # Tilføy alle deler til hovedboxen
        innhold.add(stabel_box)
        innhold.add(toga.Divider())
        innhold.add(toga.Divider())
        innhold.add(input_box)
        innhold.add(toga.Divider())
        innhold.add(resultat_box)
        self.main_window.show()

    async def hjelp(self, widget):
        hjelpedialog = toga.StackTraceDialog('Hjelpevindu', hjelpetittel, hjelpetekst)
        await self.main_window.dialog(hjelpedialog)


    async def enheter(self, widget):
        enheter = ', '.join(u.__all__[62:358])
        enheter += 'min, \', ", a, d, grad, in og Å'
        hjelpedialog = toga.StackTraceDialog('Hjelpevindu', 'Enheter du kan bruke', enheter)
        await self.main_window.dialog(hjelpedialog)


def main():
    return SymKalkulator()
