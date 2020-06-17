# Shamir Secret Sharing

A Python implementation of the Shamir Secret Sharing scheme.

Splits a given key into _w_ shares, so that at least _t_ shares are required to successfully reconstruct the key.

## Sample Usage

Create a Shamir Secret Sharing object with the chosen parameters, for example:

```
sss = SSS(13, 17, 5, 3)
```

In this example, `13` is the key, `17` is the modulo (it must be prime), `5` is the number of shares to split the key into, and `3` is the minimum number of shares required to reconstruct the key.

Generate the _x_ values, the first part of the coordinates:

```
sss.choose_x()
```

Now generate the _y_ values (the actual shares):

```
sss.generate_shares()
```

You can now reconstruct the key by passing the _x_ and _y_ coordinates to the `reconstruct_key` method:

```
print("The key is " + str(sss.reconstruct_key([1, 3, 5], [8, 10, 11])))
```

This example will return `13`, which is our key.