Программа находится в папке Source.cpp
В файл input.txt подаются на вход две строчки.
Пример:
0 0 0 2 2 0 - координаты начала и конца первого отрезка
0 0 0 2 2 0 - координаты начала и конца второго отрезка
В файл output.txt выводится результат.
Результат может быть следующим:
1. Проверка на то, что отрезки не вырождаются в точку. Если один вырождается в точку, то выводится в файл сообщение 
   об этом.
2. Проверка на то, что прямые, на которых лежат отрезки не скрещиваются. Если скрещиваются, то выводится в файл сообщение об
   этом. Если прямые сркещиваются, то они не пересекаются, и происходит выход из программы. 
3. Проверка на то, что направляющие векторы коллинеарны. Если направляющие векторы не коллинеарны, то прямые пересекаются.
   Ищется точка пересечения прямых, на которых лежат отрезки. Если точка принадлежит двум отрезкам, выводится результат. 
   Если точка не принадлежит хотя бы одному отрезку выводится информация о том, что отрезки не пересекаются
4. Если направляющие векторы коллинеарны, может быть два случая: либо они совпадают, либо параллельны. Если прямые параллельны, то отрезки не
   пересекаются, выводится сообщение об этом. Если прямые совпадают, проверяется наложение отрезков друг на друга. 