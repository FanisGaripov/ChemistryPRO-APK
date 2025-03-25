Это продолжение веб-сайта ChemistryPRO - нативные программы для ПК и Android. Все написано на pywebview. Также есть возможность прямо из под приложения перейти на веб-сайт, т.к Webview это аналог того же браузера, что позволило не переделывать сильно приложение. Компиляция(сборка) для ПК происходит в pyinstaller: pyinstaller --onefile --noconsole --add-data "templates:templates" --add-data "templates/static:templates/static" --icon=icon.ico main_for_pc.py
Компиляция для Андроид происходила в buildozer(Это было одно мучение :) ). Советую этим не заниматься самостоятельно, а скачать готовые варианты:
1) Для ПК: https://disk.yandex.ru/d/j0iTOoJV6wliAA
2) Для Android: https://disk.yandex.ru/d/7bY_p83vr9Lfbw
