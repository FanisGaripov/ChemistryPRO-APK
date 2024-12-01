import re
from chempy import balance_stoichiometry
import os
import json
import requests
from bs4 import BeautifulSoup
from urllib.parse import quote
import random
from datetime import datetime
import webview
# импортируем все библиотеки

h = ['']
pravilno = 0
otvety = 0


class ChemistryCalculator:
    def molecular_mass(self, formula):
        # Словарь с атомными массами элементов
        atomic_masses = {
            'H': 1.008,
            'He': 4.0026,
            'Li': 6.94,
            'Be': 9.0122,
            'B': 10.81,
            'C': 12.011,
            'N': 14.007,
            'O': 15.999,
            'F': 18.998,
            'Ne': 20.180,
            'Na': 22.99,
            'Mg': 24.305,
            'Al': 26.982,
            'Si': 28.085,
            'P': 30.974,
            'S': 32.06,
            'Cl': 35.45,
            'Ar': 39.948,
            'K': 39.098,
            'Ca': 40.078,
            'Sc': 44.956,
            'Ti': 47.867,
            'V': 50.941,
            'Cr': 51.996,
            'Mn': 54.938,
            'Fe': 55.845,
            'Co': 58.933,
            'Ni': 58.693,
            'Cu': 63.546,
            'Zn': 65.38,
            'Ga': 69.723,
            'Ge': 72.630,
            'As': 74.922,
            'Se': 78.971,
            'Br': 79.904,
            'Kr': 83.798,
            'Rb': 85.468,
            'Sr': 87.62,
            'Y': 88.906,
            'Zr': 91.224,
            'Nb': 92.906,
            'Mo': 95.95,
            'Tc': 98,
            'Ru': 101.07,
            'Rh': 102.905,
            'Pd': 106.42,
            'Ag': 107.868,
            'Cd': 112.414,
            'In': 114.818,
            'Sn': 118.710,
            'Sb': 121.760,
            'Te': 127.60,
            'I': 126.904,
            'Xe': 131.293,
            'Cs': 132.905,
            'Ba': 137.327,
            'La': 138.905,
            'Ce': 140.116,
            'Pr': 140.907,
            'Nd': 144.242,
            'Pm': 145,
            'Sm': 150.36,
            'Eu': 151.964,
            'Gd': 157.25,
            'Tb': 158.925,
            'Dy': 162.500,
            'Ho': 164.930,
            'Er': 167.259,
            'Tm': 168.934,
            'Yb': 173.04,
            'Lu': 174.966,
            'Hf': 178.49,
            'Ta': 180.947,
            'W': 183.84,
            'Re': 186.207,
            'Os': 190.23,
            'Ir': 192.217,
            'Pt': 195.084,
            'Au': 196.967,
            'Hg': 200.592,
            'Tl': 204.38,
            'Pb': 207.2,
            'Bi': 208.980,
            'Po': 209,
            'At': 210,
            'Rn': 222,
            'Fr': 223,
            'Ra': 226,
            'Ac': 227,
            'Th': 232.038,
            'Pa': 231.035,
            'U': 238.028,
            'Np': 237,
            'Pu': 244,
            'Am': 243,
            'Cm': 247,
            'Bk': 247,
            'Cf': 251,
            'Es': 252,
            'Fm': 257,
            'Md': 258,
            'No': 259,
            'Lr': 262,
            'Rf': 267,
            'Db': 270,
            'Sg': 271,
            'Bh': 270,
            'Hs': 277,
            'Mt': 276,
            'Ds': 281,
            'Rg': 282,
            'Cn': 285,
            'Nh': 286,
            'Fl': 289,
            'Mc': 290,
            'Lv': 293,
            'Ts': 294,
            'Og': 294,
        }

        def parse_formula(formula):
            stack = []
            current = {}
            i = 0
            while i < len(formula):
                if formula[i] == '(' or formula[i] == '[':
                    stack.append(current)
                    current = {}
                    i += 1
                elif formula[i] == ')' or formula[i] == ']':
                    i += 1
                    num = ''
                    while i < len(formula) and formula[i].isdigit():
                        num += formula[i]
                        i += 1
                    multiplier = int(num) if num else 1
                    for element, count in current.items():
                        current[element] = count * multiplier
                    if stack:
                        parent = stack.pop()
                        for element, count in current.items():
                            if element in parent:
                                parent[element] += count
                            else:
                                parent[element] = count
                        current = parent
                else:
                    match = re.match(r'([A-Z][a-z]?)(\d*)', formula[i:])
                    if match:
                        element, count = match.groups()
                        count = int(count) if count else 1
                        if element in current:
                            current[element] += count
                        else:
                            current[element] = count
                        i += len(match.group(0))
                    else:
                        i += 1
            return current

        element_counts = parse_formula(formula)
        mass = 0.0
        element_details = []
        for element, count in element_counts.items():
            element_mass = atomic_masses[element]
            total_mass = element_mass * count
            element_details.append((element, atomic_masses[element], count, total_mass))
            mass += total_mass

        return round(mass, 2), element_details


    def electronic_configuration(self, element):
        elements_data = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
            'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
            'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
            'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
            'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
            'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
            'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
            'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
            'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
            'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
            'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
            'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
            'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116,
            'Ts': 117, 'Og': 118
        }

        atomic_number = elements_data.get(element)
        atomic_number1 = elements_data.get(element)
        atom = atomic_number
        if element == '':
            return "Введите элемент", ""
        if atomic_number is None:
            return "Элемент не найден", ""

        configurations = []
        configurations1 = []
        subshells = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p', '7s', '5f',
                     '6d', '7p']
        electrons = [2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6]

        subshells1 = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f', '5s', '5p', '5d', '5f', '6s', '6p', '6d',
                      '7s', '7p']
        electrons1 = [2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 2, 6]

        for i in range(len(subshells)):
            if atomic_number > 0:
                if atomic_number >= electrons[i]:
                    configurations.append(f"{subshells[i]}^{electrons[i]}")
                    atomic_number -= electrons[i]
                else:
                    configurations.append(f"{subshells[i]}^{atomic_number}")
                    break

        for i in range(len(subshells1)):
            if atomic_number1 > 0:
                if atomic_number1 >= electrons1[i]:
                    configurations1.append(f"{subshells1[i]}^{electrons1[i]}")
                    atomic_number1 -= electrons1[i]
                else:
                    configurations1.append(f"{subshells1[i]}^{atomic_number1}")
                    break

        configuration_string = ' '.join(configurations)

        configuration_string2 = ' '.join(configurations1)

        def generate_graphical_representation(configurations):
            # графическое представление электронной конфигурации
            representation = []
            grouped_representation = {}

            for config in configurations:
                subshell, count = config.split('^')
                count = int(count)

                if subshell[0] not in grouped_representation:
                    grouped_representation[subshell[0]] = []

                cells = []

                if subshell.endswith('s'):
                    for _ in range(1):
                        if count > 0:
                            cells.append('[↑]')
                            count -= 1
                        if count > 0:
                            cells[0] += '[↓]'
                            count -= 1
                        else:
                            cells.append('[ ]')  # Пустая ячейка

                elif subshell.endswith('p'):
                    # 3 p-орбитали
                    for i in range(3):
                        if count > 0:
                            cells.append('[↑]')
                            count -= 1
                        else:
                            cells.append('[ ]')  # Пустая ячейка

                    for i in range(3):
                        if count > 0:
                            cells[i] += '[↓]'
                            count -= 1

                elif subshell.endswith('d'):
                    # 5 d-орбиталей
                    for i in range(5):
                        if count > 0:
                            cells.append('[↑]')
                            count -= 1
                        else:
                            cells.append('[ ]')  # Пустая ячейка

                    for i in range(5):
                        if count > 0:
                            cells[i] += '[↓]'
                            count -= 1

                elif subshell.endswith('f'):
                    # 7 f-орбиталей
                    for i in range(7):
                        if count > 0:
                            cells.append('[↑]')
                            count -= 1
                        else:
                            cells.append('[ ]')  # Пустая ячейка

                    for i in range(7):
                        if count > 0:
                            cells[i] += '[↓]'
                            count -= 1

                grouped_representation[subshell[0]].append(f"{subshell}: " + ' '.join(cells))

            # Сборка финального представления
            for level in sorted(grouped_representation.keys()):
                representation.extend(grouped_representation[level])

            return "\n".join(representation)

        # Графическое представление(текстовое, используются [↑] и [↓], открывающая и закрывающая скобка - это одна клетка)
        graphic_representation = generate_graphical_representation(configurations)

        return {
            "configuration_string": configuration_string2,
            "configuration_string2": configuration_string2,
            "graphic_representation": graphic_representation,
            "atom": elements_data[element]
        }


    def uravnivanie(self, formula):
        # баланс уравнений
        reactants_input, products_input = formula.split('=')
        reactants = {x.split()[0].strip(): int(x.split()[1]) if len(x.split()) > 1 else 1 for x in
                     reactants_input.split('+')}
        products = {x.split()[0].strip(): int(x.split()[1]) if len(x.split()) > 1 else 1 for x in products_input.split('+')}

        balanced_reaction = balance_stoichiometry(reactants, products)

        reactants_str = ' + '.join([f"{v}{k}" if v != 1 else f"{k}" for k, v in balanced_reaction[0].items()])
        products_str = ' + '.join([f"{v}{k}" if v != 1 else f"{k}" for k, v in balanced_reaction[1].items()])

        otvet = f"{reactants_str} = {products_str}"

        return otvet


    def get_chemical_equation_solution(self, reaction):
        '''метод обработчик дописывания хим.реакций. коротко о нем: принимает из основной функции реакцию, вставляет
        ее в ссылку и возвращает ответ, который парсит(выкидывает все лишнее) только до нужных строчек'''
        # Кодируем реакцию для URL
        encoded_reaction = quote(reaction)

        # Формируем URL с учетом химической реакции
        url = f"https://chemequations.com/ru/?s={encoded_reaction}"

        # Отправляем GET-запрос
        response = requests.get(url)

        # Проверка успешности запроса
        if response.status_code == 200:
            # Парсим HTML-ответ
            soup = BeautifulSoup(response.text, 'html.parser')

                # Находим элемент с классом "equation main-equation well"
            result = soup.find('h1', class_='equation main-equation well')

            if result:
                react1 = result.get_text(strip=True)
                if '(g)' in react1:
                    react1 = react1.replace('(g)', '')
                if '(s)' in react1:
                    react1 = react1.replace('(s)', '')
                if '(aq)' in react1:
                    react1 = react1.replace('(aq)', '')
                if '(l)' in react1:
                    react1 = react1.replace('(l)', '')
                return react1
            else:
                return 'Решение не найдено.'
        else:
            return f"Ошибка при запросе: {response.status_code}"


    def get_reaction_chain(self, reaction):
        # цепочка превращений
        url = f"https://chemer.ru/services/reactions/chains/{reaction}"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
        }
        session = requests.Session()
        session.headers.update(headers)
        response = session.get(url)

        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            inset_divs = soup.find_all('div', class_='inset')  # Ищем все группы
            reaction_groups = {}  # Словарь для групп реакций

                # Сначала собираем группы из h2
            for inset_div in inset_divs:
                product_header = inset_div.find('h2')
                if product_header:
                    product_text = product_header.get_text().strip()
                    reaction_groups[product_text] = []  # Создаем пустой список для реакций этой группы

            # Теперь собираем реакции и добавляем их в соответствующие группы
            content_sections = soup.find_all('section', class_='content')  # Ищем все секции с классом 'content'
            if content_sections:
                for content_section in content_sections:
                    reactions = content_section.find_all('p', class_='resizable-block')  # Ищем все 'p' внутри каждой секции
                    if reactions:
                        # Находим соответствующую группу для текущей секции
                        group_header = content_section.find_previous('div', class_='inset').find('h2')
                        if group_header:
                            group_name = group_header.get_text().strip()
                            for reaction in reactions:
                                reaction_text = reaction.get_text().strip()  # Извлекаем текст реакции
                                if group_name in reaction_groups:
                                    reaction_groups[group_name].append(reaction_text)  # Добавляем реакцию в соответствующую группу

                # Формируем окончательный результат
            final_results = []
            for group, reactions in reaction_groups.items():
                final_results.append(f"Как из {group}:")
                final_results.extend(reactions)  # Добавляем все реакции для данной группы

            return final_results  # Возвращаем сгруппированные реакции
        else:
            return [f"Ошибка при запросе: {response.status_code}"]  # Обработка ошибки запроса
        return []



    def rastvory(self, mass_solution, mass_solute, mass_fraction):
        # калькулятор растворимостей
        # Выполняем расчеты
        if mass_solution and mass_solute is not None:
            # Если известны масса раствора и масса вещества, рассчитываем массовую долю
            mass_fraction = (mass_solute / mass_solution) * 100
    
        elif mass_solution and mass_fraction is not None:
            # Если известны масса раствора и массовая доля, рассчитываем массу растворенного вещества
            mass_solute = (mass_fraction / 100) * mass_solution
    
        elif mass_solute and mass_fraction is not None:
            # Если известны масса вещества и массовая доля, рассчитываем массу раствора
            mass_solution = mass_solute / (mass_fraction / 100)
    
        return mass_solution, mass_solute, mass_fraction


    def minigame(self, element):
        def minigamefunc():
            # функция обработчик миниигры
            a = random.randint(0, 117)
            atomic_masses = {
                'H': 'Водород',
                'He': 'Гелий',
                'Li': 'Литий',
                'Be': 'Бериллий',
                'B': 'Бор',
                'C': 'Углерод',
                'N': 'Азот',
                'O': 'Кислород',
                'F': 'Фтор',
                'Ne': 'Неон',
                'Na': 'Натрий',
                'Mg': 'Магний',
                'Al': 'Алюминий',
                'Si': 'Кремний',
                'P': 'Фосфор',
                'S': 'Сера',
                'Cl': 'Хлор',
                'Ar': 'Аргон',
                'K': 'Калий',
                'Ca': 'Кальций',
                'Sc': 'Скандий',
                'Ti': 'Титан',
                'V': 'Ванадий',
                'Cr': 'Хром',
                'Mn': 'Марганец',
                'Fe': 'Железо',
                'Co': 'Кобальт',
                'Ni': 'Никель',
                'Cu': 'Медь',
                'Zn': 'Цинк',
                'Ga': 'Галлий',
                'Ge': 'Германий',
                'As': 'Мышьяк',
                'Se': 'Селен',
                'Br': 'Бром',
                'Kr': 'Криптон',
                'Rb': 'Рубидий',
                'Sr': 'Стронций',
                'Y': 'Иттрий',
                'Zr': 'Цирконий',
                'Nb': 'Ниобий',
                'Mo': 'Молибден',
                'Tc': 'Технеций',
                'Ru': 'Рутений',
                'Rh': 'Родий',
                'Pd': 'Палладий',
                'Ag': 'Серебро',
                'Cd': 'Кадмий',
                'In': 'Индий',
                'Sn': 'Олово',
                'Sb': 'Сурьма',
                'Te': 'Теллур',
                'I': 'Йод',
                'Xe': 'Ксенон',
                'Cs': 'Цезий',
                'Ba': 'Барий',
                'La': 'Лантан',
                'Ce': 'Церий',
                'Pr': 'Празеодим',
                'Nd': 'Неодим',
                'Pm': 'Прометий',
                'Sm': 'Самарий',
                'Eu': 'Европий',
                'Gd': 'Гадолиний',
                'Tb': 'Тербий',
                'Dy': 'Диспрозий',
                'Ho': 'Гольмий',
                'Er': 'Эрбий',
                'Tm': 'Тулий',
                'Yb': 'Иттербий',
                'Lu': 'Лютеций',
                'Hf': 'Гафний',
                'Ta': 'Тантал',
                'W': 'Вольфрам',
                'Re': 'Рений',
                'Os': 'Осмий',
                'Ir': 'Иридий',
                'Pt': 'Платина',
                'Au': 'Золото',
                'Hg': 'Ртуть',
                'Tl': 'Таллий',
                'Pb': 'Свинец',
                'Bi': 'Висмут',
                'Po': 'Полоний',
                'At': 'Астат',
                'Rn': 'Радон',
                'Fr': 'Франций',
                'Ra': 'Радий',
                'Ac': 'Актиний',
                'Th': 'Торий',
                'Pa': 'Проактиний',
                'U': 'Уран',
                'Np': 'Нептуний',
                'Pu': 'Плутоний',
                'Am': 'Америций',
                'Cm': 'Кюрий',
                'Bk': 'Берклий',
                'Cf': 'Калифорний',
                'Es': 'Эйнштейний',
                'Fm': 'Фермий',
                'Md': 'Менделевий',
                'No': 'Нобелий',
                'Lr': 'Лоуренсий',
                'Rf': 'Резерфордий',
                'Db': 'Дубний',
                'Sg': 'Сиборгий',
                'Bh': 'Борий',
                'Hs': 'Хассий',
                'Mt': 'Майтнерий',
                'Ds': 'Дармштадтий',
                'Rg': 'Рентгений',
                'Cn': 'Коперниций',
                'Nh': 'Нихоний',
                'Fl': 'Флеровий',
                'Mc': 'Московий',
                'Lv': 'Ливерморий',
                'Ts': 'Теннессин',
                'Og': 'Оганессон',
            }
            k = []
            d = ''
            b = ""
            for i in atomic_masses:
                k.append(i)
            b = k[a]
            nazv = atomic_masses[b]
            print(b, nazv)
            return b, nazv
        '''функция, которая возвращает страницу мини-игры
        Коротко о мини-игре:
        Это Игра для запоминания элементов таблицы Менделеева.
        Выводится элемент, а игрок должен написать, то как он называется на РУССКОМ языке'''
        d = ""
        global h, pravilno, otvety
        res = minigamefunc()
        b = res[0]
        nazv = res[1]
        h.append(nazv)
        if res:
            if element == h[-2]:
                d = 'Верно, следующий'
                pravilno += 1
                otvety += 1
                if pravilno == 10:
                    right_percent = round((pravilno / otvety) * 100, 2)
                    pravilno = 0
                    otvety = 0
                    if right_percent >= 50:
                        print('Вы победили')
                    else:
                        print('Вы проиграли')
                    h = ['']
            else:
                d = f'Неправильно, ответ: {h[-2]}'
                otvety += 1
            return {
                "b": b,
                "d": d,
                "pravilno": pravilno,
                "otvety": otvety
            }


    def get_substance_html(self, substance_name):
        # получение имени орг вещества из таблицы на сайте при помощи парсинга этой страницы
        url = "https://chemer.ru/services/organic/structural"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
        }
        session = requests.Session()
        session.headers.update(headers)
        response = session.get(url)
        global klass
        klass = ''
        namez = []

        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            table = soup.find('table')
            rows = table.find_all('tr')

            for row in rows:
                cols = row.find_all('td')
                if cols:
                    name = cols[0].text.strip()
                    klass = cols[1].text.strip()
                    link = cols[0].find('a')['href']
                    if substance_name.lower() in name.lower() and substance_name.lower() == name.lower():
                        substance_url = f"https://chemer.ru/services/organic/{link}"
                        substance_response = session.get(substance_url)
                        return substance_response.text, None
                    elif substance_name.lower() in name.lower() and name.lower()[2:] != substance_name.lower():
                        namez.append(name)
                    elif substance_name.lower() in name.lower() and name.lower()[:2] == 'н-' and name.lower()[2:] == substance_name.lower():
                        substance_url = f"https://chemer.ru/services/organic/{link}"
                        substance_response = session.get(substance_url)
                        return substance_response.text, None
        return None, namez


    def extract_svg_and_symbols(self, html_code):
        soup = BeautifulSoup(html_code, 'html.parser')
        svg_elements = soup.find_all('svg')
        symbols = soup.find_all('symbol')
        names = ''

        if not svg_elements:
            return None, None, None

        first_svg_content = str(svg_elements[0])
        if 'width' not in first_svg_content or 'height' not in first_svg_content:
            first_svg_content = first_svg_content.replace('<svg', '<svg width="200" height="200"')

        isomer_svgs = []
        spacing = 220  # Расстояние между изомерами
        max_per_row = 20  # Максимум изомеров в строке

        # Извлечение секции с id='tab1'
        tab1_section = soup.find('section', id='tab1')

        if tab1_section:
            svg_elements2 = tab1_section.find_all('svg')
            names = tab1_section.find_all('a')

            for index, svg in enumerate(svg_elements2):
                row = index // max_per_row
                col = index % max_per_row
                x = col - 1  # Установка x координаты
                y = row  # Установка y координаты
                svg_str = str(svg).replace('<svg', f'<svg x="{x}" y="{y}"')
                isomer_svgs.append(svg_str)

        isomer_svgs_content = ''.join(isomer_svgs)
        symbol_content = ''.join(str(symbol) for symbol in symbols)

        return first_svg_content, isomer_svgs_content, symbol_content, names

    def orghim(self, substance_name):
        def get_substance_html(substance_name):
            url = "https://chemer.ru/services/organic/structural"
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
            }
            session = requests.Session()
            session.headers.update(headers)
            response = session.get(url)
            namez = []

            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')
                table = soup.find('table')
                rows = table.find_all('tr')

                for row in rows:
                    cols = row.find_all('td')
                    if cols:
                        global klass
                        name = cols[0].text.strip()
                        klass = cols[1].text.strip()
                        link = cols[0].find('a')['href']
                        if substance_name.lower() == name.lower():
                            substance_url = f"https://chemer.ru/services/organic/{link}"
                            substance_response = session.get(substance_url)
                            return substance_response.text, None
                        elif substance_name.lower() in name.lower():
                            namez.append(name)
            return None, namez

        html_code, variants = get_substance_html(substance_name)

        if html_code:
            first_svg, isomers_svg, symbols_svg, names = self.extract_svg_and_symbols(html_code)

            # Сохраняем первую SVG-картинку и символы в файл
            if first_svg:
                svg_file_path = 'static/output.svg'
                with open('templates/static/output.svg', 'w', encoding='utf-8') as f:
                    f.write(
                        f"<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>{symbols_svg}{first_svg}</svg>")

            # Сохраняем изомеры в отдельные файлы
            isomer_files = []
            combined = []
            if isomers_svg.strip():
                for index, svg in enumerate(isomers_svg.split('</svg>')):
                    if index <= 9 and svg.strip():
                        file_name = f'templates/static/isomer_{index}.svg'
                        with open(file_name, 'w', encoding='utf-8') as f:
                            f.write(f"<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>{symbols_svg}{svg}</svg></svg>")
                        isomer_files.append(f'isomer_{index}.svg')

                for i, nazv in enumerate(names):
                    if i <= 9:
                        soup = BeautifulSoup(str(nazv), 'html.parser')
                        nazv = soup.a.text
                        combined.append((nazv.capitalize(), isomer_files[i]))

            return {
                "svg_file": 'static/output.svg',
                "isomer_files": isomer_files,
                "substance_name": substance_name,
                "klass": klass,
                "combined": combined,
                "variants": variants
            }

        return {
            "svg_file": None,
            "isomer_files": [],
            "substance_name": substance_name,
            "klass": None,
            "combined": [],
            "variants": variants
        }


# Функция для возврата на предыдущую страницу
def go_back(window):
    window.evaluate('window.history.back()')


# Функция для перехода на следующую страницу
def go_forward(window):
    window.evaluate('window.history.forward()')


def create_window():
    icon_path = 'static/favicon.svg'
    calculator = ChemistryCalculator()
    window = webview.create_window('ChemistryPRO-APP', "templates/main.html", js_api=calculator)
    webview.start()

if __name__ == '__main__':
    create_window()