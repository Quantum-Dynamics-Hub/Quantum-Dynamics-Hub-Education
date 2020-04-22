---
title: "Algorithms and methodology prototyping. Introducing Libra"
teaching: 60
exercises: 10
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

> ## Warning
>
> This webpage is still under development. Please check for updates.
>
{:.discussion}


# Tutorial 1: Deep vs. Shallow copying in Python

Perhaps, one of the most "dangerous" things you can do in Python is to copy variables.

Most of the time, copying in Python is made by reference, which is also called *shallow* copying, and this is something is so simple to overlook or disregard.

This type of copying applies to the Python objects, but not to elementary data types, such as floats or so.

So, when you are doing the following:


> x = 10
> y = x
> print(F"x = {x}, y = {y}")

> x = 20
> print(F"x = {x}, y = {y}")


{% include links.md %}

