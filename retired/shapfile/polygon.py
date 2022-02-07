from shapely.ops import polygonize, polygonize_full
lines = [
    ((0, 0),(0.5, 0.5), (1, 1),(1.5,0.5)),
    ((0, 0), (0, -1),(1, 0), (1.5,0.5)),
    ((0.5, 0.5), (1, 0))
    ]


print(lines)
print(type(lines))
print(type(lines[0]))
a= lines[0]
b=a[0]

c = polygonize(lines)
d = list(c)
print(d)