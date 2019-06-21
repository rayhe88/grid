"""Onedgrid tests file."""
from unittest import TestCase

from grid.onedgrid import (
    GaussChebyshev,
    GaussChebyshevLobatto,
    GaussChebyshevType2,
    GaussLaguerre,
    GaussLegendre,
    RectangleRuleSine,
    RectangleRuleSineEndPoints,
    TanhSinh,
    generate_onedgrid,
)

import numpy as np

from scipy.special import roots_chebyu, roots_genlaguerre, roots_legendre


class TestOneDGrid(TestCase):
    """OneDGrid test class."""

    def setUp(self):
        """Test setup function."""
        ...

    def test_generate_ondgrid(self):
        """Place holder tests for generate_onedgrid."""
        generate_onedgrid(0)

    def test_gausslaguerre(self):
        """Test Guass Laguerre polynomial grid."""
        points, weights = np.polynomial.laguerre.laggauss(10)
        roots_genlaguerre(10, 0)
        grid = GaussLaguerre(10)
        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_gausslengendre(self):
        """Test Guass Lengendre polynomial grid."""
        points, weights = roots_legendre(10)
        grid = GaussLegendre(10)
        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_gausschebyshev(self):
        """Test Guass Chebyshev polynomial grid."""
        points, weights = np.polynomial.chebyshev.chebgauss(10)
        grid = GaussChebyshev(10)
        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_gausschebyshev2(self):
        """Test Gauss Chebyshev type 2 polynomial grid."""
        points, weights = roots_chebyu(10)
        grid = GaussChebyshevType2(10)
        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_gausschebyshevlobatto(self):
        """Test Gauss Chebyshev Lobatto grid."""
        grid = GaussChebyshevLobatto(10)

        idx = np.arange(10)
        weights = np.ones(10)
        idx = (idx * np.pi) / 9

        points = np.cos(idx)
        points = np.sort(points)

        weights *= np.pi / 9
        weights[0] /= 2
        weights[9] /= 2

        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_rectanglesineendpoints(self):
        """Test for rectangle rule for sine series with endpoints."""
        grid = RectangleRuleSineEndPoints(10)

        idx = np.arange(10) + 1
        points = idx / 11

        weights = np.zeros(10)

        index_m = np.arange(10) + 1

        for i in range(0, 10):
            elements = np.zeros(10)
            elements = np.sin(index_m * np.pi * points[i])
            elements *= (1 - np.cos(index_m * np.pi)) / (index_m * np.pi)

            weights[i] = (2 / (11)) * np.sum(elements)
        
        points = 2 * points - 1

        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_rectanglesine(self):
        """Test for rectangle rule for sine series without endpoint."""
        grid = RectangleRuleSine(10)

        idx = np.arange(10) + 1
        points = (2 * idx - 1) / 20

        weights = np.zeros(10)

        index_m = np.arange(9) + 1

        for i in range(0, 10):
            elements = np.zeros(9)
            elements = np.sin(index_m * np.pi * points[i])
            elements *= np.sin(index_m * np.pi / 2) ** 2
            elements /= index_m

            weights[i] = (4 / (10 * np.pi)) * np.sum(elements)

            weights[i] += (
                (2 / (10 * np.pi ** 2))
                * np.sin(10 * np.pi * points[i])
                * np.sin(10 * np.pi / 2) ** 2
            )

        points = 2 * points - 1

        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_tanhsinh(self):
        """Test for Tanh - Sinh rule."""
        delta = 0.1 * np.pi / np.sqrt(11)
        grid = TanhSinh(11, delta)

        jmin = -5
        points = np.zeros(11)
        weights = np.zeros(11)

        for i in range(0, 11):
            j = jmin + i
            arg = np.pi * np.sinh(j * delta) / 2

            points[i] = np.tanh(arg)

            weights[i] = np.pi * delta * np.cosh(j * delta) * 0.5
            weights[i] /= np.cosh(arg) ** 2

        assert np.allclose(grid.points, points)
        assert np.allclose(grid.weights, weights)

    def test_errors_raise(self):
        """Test errors raise."""
        with self.assertRaises(ValueError):
            GaussLaguerre(10, -1)

        with self.assertRaises(ValueError):
            TanhSinh(10, 1)
