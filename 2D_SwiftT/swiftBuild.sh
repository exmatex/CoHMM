#!/bin/bash
tclsh make-package.tcl > pkgIndex.tcl
stc -r ./ 2D_SwiftT.swift
